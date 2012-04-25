package ca.mcgill.mcb.pcingola.osCmd;

import java.io.File;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Run an OS command 
 * 
 * @author pcingola
 */
public class OsCmdRunner {

	boolean executing = false, started = false;
	int exitValue = 0;
	long defaultWaitTime = 100; // Default time to use in 'wait' calls when initialting the command
	long defaultLoopWaitTime = 1000; // Default time to use in 'wait' calls
	String jobId = "";
	String head = ""; // Ouput's head
	String headStderr = ""; // Stderr's head
	String error; // Latest error message
	ExecuteOsCommand osCmd = null;

	/**
	* Run an OS command only if the output files does not exists.
	* 
	* 	opts[0] 			: OS Command 
	* 	opts[1] ... opts[N]	: Command line options
	* 	outputFile			: Where the results are stored (if the file exists, the command is NOT run)
	* 	redirect			: If 'redirect=true' then run "command > outputFile" (i.e. redirect STDOUT to 'outputFile'). Output is assumed to be binary.
	* 
	* @param opts
	* @param outputFile
	* @param redirectToOutput
	* @return true if command excecuted OK or outputFile exists
	*/
	public static boolean runIfNotExists(String[] opts, String outputFile, boolean redirectToOutput) {
		// Already done?
		if (Gpr.exists(outputFile)) return true;

		try {
			// Create command and execute it
			String id = opts[0];
			OsCmdRunner cmd = new OsCmdRunner(id, opts);
			cmd.getOsCmd().setQuiet(false, false);

			// Have to redirect?
			if (redirectToOutput) {
				cmd.getOsCmd().setBinaryStdout(true);
				cmd.getOsCmd().setRedirectStdout(outputFile);
			}

			Timer.show("\tExecuting command: " + cmd.getOsCmd());
			cmd.run();
		} catch (Throwable t) {
			t.printStackTrace();

			// Try to delete output file
			if ((outputFile != null) && (!outputFile.isEmpty())) {
				try {
					(new File(outputFile)).delete();
				} catch (Exception e) {
					// Nothing to do
				}
			}

			return false;
		}

		return true;
	}

	public OsCmdRunner(String jobId, String osCmdStr[]) {
		super();
		this.jobId = jobId;
		osCmd = new ExecuteOsCommand(osCmdStr);
	}

	/**
	 * Close (kill) command
	 */
	public synchronized void close() {
		if (osCmd != null) {
			osCmd.kill();
			head = osCmd.getHead();
			headStderr = osCmd.getHeadStderr();
		}
		osCmd = null;
	}

	/**
	 * Stop execution of this thread
	 */
	public synchronized void finish() {
		executing = false; // Set run to false and wake up from 'wait'. See run() method
		notify();
	}

	public long getDefaultWaitTime() {
		return defaultWaitTime;
	}

	public String getError() {
		return error;
	}

	public int getExitValue() {
		return exitValue;
	}

	public String getHead() {
		return head;
	}

	public String getHeadStderr() {
		return headStderr;
	}

	public String getJobId() {
		return jobId;
	}

	public ExecuteOsCommand getOsCmd() {
		return osCmd;
	}

	public int getProgress() {
		return (osCmd == null ? 0 : osCmd.getProgress());
	}

	/**
	 * Has this runner finished?
	 * @return
	 */
	public boolean isDone() {
		return started && !executing; // In order to finish, it has to be started and not running any more
	}

	public boolean isExecuting() {
		if (osCmd == null) return false;
		return executing && osCmd.isExecuting();
	}

	public void run() {
		try {
			executing = true; // OK, we are running
			osCmd.start();

			Thread.sleep(defaultWaitTime); // Allow some time for the thread to start

			// Wait for command to start 
			while (!osCmd.isStarted() && isExecuting())
				Thread.sleep(defaultWaitTime);

			// Wait for stdin to became available
			while ((osCmd.getStdin() == null) && isExecuting())
				Thread.sleep(defaultWaitTime);

			// Now the comamnd started executing
			started = true;

			synchronized (this) {
				osCmd.setObjetcToNotify(this); // Notify me when done (i.e. the command finished)
				while (isExecuting())
					wait(defaultLoopWaitTime);
			}

		} catch (Throwable t) {
			error = t.toString();
			t.printStackTrace(); // Something happened? => Stop this thread
		} finally {
			exitValue = osCmd.getExitValue();
			close();
			executing = false;
			started = true;
		}
	}

	public void setDefaultWaitTime(long defaultWaitTime) {
		this.defaultWaitTime = defaultWaitTime;
	}
}
