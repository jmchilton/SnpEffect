package ca.mcgill.mcb.pcingola.stats.plot;

public class GoogleBarChart extends GoogleLineChart {

	public GoogleBarChart(String name) {
		super(name);
	}

	public GoogleBarChart(String name, int width, int height) {
		super(name, width, height);
	}

	@Override
	public String toStringHtmlHeader() {
		StringBuilder sb = new StringBuilder();
		sb.append("<script type=\"text/javascript\" src=\"http://www.google.com/jsapi\"></script>");
		sb.append("<script type=\"text/javascript\"> google.load('visualization', '1', {packages: ['corechart']}); </script>\n");
		sb.append("<script type=\"text/javascript\">\n");
		sb.append("\tfunction draw_" + id + "() {\n");
		sb.append("\t\tvar data = google.visualization.arrayToDataTable([\n");

		// Column titles
		sb.append("\t[ '' , ");
		int i = 0;
		for (String ct : columnTitltes) {
			sb.append((i > 0 ? "," : "") + "'" + ct + "'");
			i++;
		}
		sb.append("]\n");

		// Data
		int maxLen = maxColumnLength();
		for (i = 0; i < maxLen; i++) {
			// X labels
			String lab = getXLabel(i);
			if (lab != null) lab = "'" + lab + "'";
			sb.append("\t,[ " + lab);

			// Data
			for (int j = 0; j < columns.size(); j++)
				sb.append("," + getValue(i, j));
			sb.append("]\n");
		}

		sb.append("\t\t]);\n");
		sb.append("\t\tvar ac = new google.visualization.ColumnChart(document.getElementById('visualization_" + id + "'));\n");
		sb.append("\t\tac.draw(data, { title : '" + title + "', isStacked: " + stacked + ", width: " + width + ", height: " + height + ", vAxis: {title: \"" + vAxis + "\"}, hAxis: {title: \"" + hAxis + "\"} });\n");
		sb.append("\t\t}\n");
		sb.append("\tgoogle.setOnLoadCallback(draw_" + id + ");\n");
		sb.append("</script>\n");
		sb.append("\n");
		return sb.toString();
	}
}