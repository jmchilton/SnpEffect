package ca.mcgill.mcb.pcingola.stats.plot;

import java.util.ArrayList;

import ca.mcgill.mcb.pcingola.util.Gpr;

public class GoogleLineChart {

	boolean stacked = false;

	int width = 800, height = 600;
	String name;
	String title = "";
	String vAxis = "";
	String hAxis = "";
	ArrayList<String> columnTitltes;
	ArrayList<String> xLables;
	ArrayList<ArrayList<String>> columns;

	public GoogleLineChart(String name) {
		this.name = name;
		init();
	}

	public GoogleLineChart(String name, int width, int height) {
		this.name = name;
		this.width = width;
		this.height = height;
		init();
	}

	public void addColumn(String colTitle, ArrayList<String> columnValues) {
		columns.add(columnValues);
		columnTitltes.add(colTitle);
	}

	String getValue(int i, int j) {
		if (j >= columns.size()) return null;
		ArrayList<String> col = columns.get(j);
		if (i >= col.size()) return null;
		return col.get(i);
	}

	String getXLabel(int idx) {
		if (idx >= xLables.size()) return null;
		return xLables.get(idx);
	}

	void init() {
		columnTitltes = new ArrayList<String>();
		xLables = new ArrayList<String>();
		columns = new ArrayList<ArrayList<String>>();
	}

	/**
	 * Max column length
	 * @return
	 */
	int maxColumnLength() {
		int size = 0;
		for (int i = 0; i < columns.size(); i++) {
			size = Math.max(size, columns.get(i).size());
			Gpr.debug("\t\tColumn : " + i + "\tlength: " + columns.get(i).size());

		}
		Gpr.debug("Max col length: " + size);
		return size;
	}

	public void sethAxis(String hAxis) {
		this.hAxis = hAxis;
	}

	public void setHeight(int height) {
		this.height = height;
	}

	public void setName(String name) {
		this.name = name;
	}

	public void setStacked(boolean stacked) {
		this.stacked = stacked;
	}

	public void setTitle(String title) {
		this.title = title;
	}

	public void setvAxis(String vAxis) {
		this.vAxis = vAxis;
	}

	public void setWidth(int width) {
		this.width = width;
	}

	public void setxLables(ArrayList<String> xLables) {
		this.xLables = xLables;
	}

	public String toStringHtmlBody() {
		StringBuilder sb = new StringBuilder();
		sb.append("<script type=\"text/javascript\" src=\"http://www.google.com/jsapi\"></script>");
		sb.append("<script type=\"text/javascript\"> google.load('visualization', '1', {packages: ['corechart']}); </script>\n");
		sb.append("<script type=\"text/javascript\">\n");
		sb.append("\tfunction draw_" + name + "() {\n");
		sb.append("\t\tvar data = google.visualization.arrayToDataTable([\n");

		// Column titles
		sb.append("\t[ '' , ");
		int i = 0;
		for (String ct : columnTitltes) {
			sb.append((i > 0 ? "," : "") + "'" + ct + "'");
			i++;
		}
		sb.append("]\n");

		// Date
		int maxLen = maxColumnLength();
		for (i = 0; i < maxLen; i++) {
			sb.append("\t,[ " + getXLabel(i));
			for (int j = 0; j < columns.size(); j++)
				sb.append("," + getValue(i, j));
			sb.append("]\n");
		}

		sb.append("\t\t]);\n");
		sb.append("\t\tvar ac = new google.visualization.AreaChart(document.getElementById('visualization_" + name + "'));\n");
		sb.append("\t\tac.draw(data, { title : '" + title + "', isStacked: " + stacked + ", width: " + width + ", height: " + height + ", vAxis: {title: \"" + vAxis + "\"}, hAxis: {title: \"" + hAxis + "\"} });\n");
		sb.append("\t\t}\n");
		sb.append("\tgoogle.setOnLoadCallback(draw_" + name + ");\n");
		sb.append("</script>\n");
		sb.append("\n");
		return sb.toString();
	}

	public String toStringHtmlHeader() {
		StringBuilder sb = new StringBuilder();
		sb.append("<div id=\"visualization_" + name + "\" style=\"width: " + width + "px; height: " + width + "px;\"></div>\n");
		return sb.toString();
	}

}
