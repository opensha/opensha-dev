package scratch.kevin;

import java.util.Arrays;

public class TracTableCreator {
	
	public static String getTracTableString(Object[][] table) {
		return getTracTableString(Arrays.asList(table));
	}
	
	public static String getTracTableString(Iterable<? extends Object[]> table) {
		String str = null;
		
		for (Object[] row : table) {
			if (str == null)
				str = "";
			else
				str += "\n";
			for (Object cell : row) {
				str += "||"+cell;
			}
			str += "||";
		}
		
		return str;
	}

}
