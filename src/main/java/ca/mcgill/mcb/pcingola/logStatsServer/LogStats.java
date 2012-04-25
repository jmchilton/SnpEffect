package ca.mcgill.mcb.pcingola.logStatsServer;

import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Properties;

import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.util.Gpr;

/** 
 * 
 * Log basic usage information to a server (for feedback and stats)
 * This information an always be suppressed (no info sent at all)
 * 
 */
public class LogStats extends Thread {

	public enum RequestResult {
		OK, ERROR, NOINFO;

		/** A request is "Completed" if the communication worked */
		public boolean completed() {
			return this != ERROR && this != NOINFO;
		}
	}

	public static boolean debug = false; // Debug mode?

	// Log server parameters
	private static final String URL_ROOT = "http://www.tacner.com/special/recuse.php";
	private static final String HTTP_CHARSET = "ISO-8859-1";
	private static final int HTTP_CONNECT_TIMEOUT_MSECS = 22000;
	private static final int HTTP_READ_TIMEOUT_MSECS = 23000;

	// Connection
	private final Properties response = new Properties(); // empty if bad connection

	// Protocol
	public StringBuilder msg = new StringBuilder(); // info for the user
	private final String version;

	private RequestResult res = RequestResult.NOINFO;
	private long duration; // time to complete the request, in msecs - succuessfull or not

	HashMap<String, String> values;

	// Constructor
	public LogStats() {
		version = SnpEff.VERSION_SHORT;
		values = new HashMap<String, String>();
	}

	/**
	 * Add a 'name=value' pair
	 * @param name
	 * @param value
	 */
	public void add(String name, String value) {
		values.put(name, value);
	}

	/** Unencripted */
	private URL buildUrl() {
		try {
			return buildUrl0();
		} catch (MalformedURLException e) {
			msg.append(e.getMessage());
			return null;
		}
	}

	/**
	 * Create URL
	 * @return
	 * @throws MalformedURLException
	 */
	private URL buildUrl0() throws MalformedURLException {
		StringBuilder urlsb = new StringBuilder();
		urlsb.append(URL_ROOT).append("?");

		// Add program and version
		urlsb.append("program=").append(encode2url(SnpEff.class.getSimpleName()));
		urlsb.append("&version=").append(encode2url(version));

		// Add all other 'name=value' pairs in alphabetical order
		ArrayList<String> names = new ArrayList<String>();
		names.addAll(values.keySet());
		Collections.sort(names);
		for (String name : names)
			urlsb.append("&" + name + "=").append(encode2url(values.get(name)));

		return new URL(urlsb.toString());
	}

	public void connect() {
		// Step 0 : Build URL
		int step = 0;
		long t0 = System.currentTimeMillis();
		if (debug) Gpr.debug("Connect Step = " + step);
		try {
			URL url = buildUrl();

			// Step 1: Open connection
			step = 1;
			if (debug) Gpr.debug("Connect Step = " + step);
			URLConnection hc = url.openConnection();

			// Step 2: Set parameters
			step = 2;
			if (debug) Gpr.debug("Connect Step = " + step);
			hc.setConnectTimeout(HTTP_CONNECT_TIMEOUT_MSECS);
			hc.setReadTimeout(HTTP_READ_TIMEOUT_MSECS);

			// Step 3: Connect to server
			step = 3;
			if (debug) Gpr.debug("Connect Step = " + step);
			response.load(hc.getInputStream());

			// Step 4: Parse results (nothing done here)
			step = 4;
			if (debug) Gpr.debug("Connect Step = " + step);
			res = RequestResult.OK;
		} catch (Exception e) {
			msg.append(step > 3 ? "Bad response" : "Error in connection. ").append(" Step " + step).append("(").append(e.toString()).append(")");
			res = RequestResult.ERROR;
		} finally {
			duration = System.currentTimeMillis() - t0;
			if (debug && !res.completed()) Gpr.debug("Error in connection: " + res + " step=" + step + " duration(msecs)=" + duration + " " + msg);
		}

		if (debug) Gpr.debug("Connect done!");
	}

	/**
	 * Encode data to URL format
	 * @param data
	 * @return
	 */
	private String encode2url(String data) {
		try {
			return URLEncoder.encode(data, HTTP_CHARSET);
		} catch (UnsupportedEncodingException e) {
			return "";
		}
	}

	public RequestResult getRes() {
		return res;
	}

	/**
	 * Run thread in background
	 */
	@Override
	public void run() {
		try {
			if (debug) Gpr.debug("Running thread");
			connect();
			if (debug) Gpr.debug("Thread finished");
		} catch (Throwable t) {
			if (debug) t.printStackTrace();; // Do nothing if it fails
		}

	}
}
