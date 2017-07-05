package scratch.kevin.cybershake.etasCalcs;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class StampedeSpeedPlotter {

	/**
	 * this will plot completion times as a function of threads for MPJTaskCalculator logs
	 * 
	 *  assumes directory name contains "-#thread" for parsing number of threads
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File("/tmp/stampede_log");
		
		SimpleDateFormat df = MPJTaskCalculator.df;
		
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		
		for (File subDir : dir.listFiles()) {
			if (!subDir.isDirectory())
				continue;
			String name = subDir.getName();
			if (name.startsWith("."))
				continue;
			name = name.substring(0, name.indexOf("thread"));
			name = name.substring(name.lastIndexOf("-")+1);
			int threads = Integer.parseInt(name);
			
			for (File file : subDir.listFiles()) {
				name = file.getName();
				if (name.contains(".pbs.o")) {
					Date firstDate = null;
					Date lastDate = null;
					
					BufferedReader tis =
							new BufferedReader(new InputStreamReader(new FileInputStream(file)));
					String str = tis.readLine();
					while(str != null) {
						if (str.startsWith("[")) {
							// isolate date
							String dateStr = str.substring(1);
							dateStr = dateStr.substring(0, dateStr.indexOf(" "));
							try {
								Date date = df.parse(dateStr);
								if (firstDate == null)
									firstDate = date;
//								if (str.contains("Dispatcher") && str.contains("DONE"))
									lastDate = date;
							} catch (ParseException e) {
								// do nothing on bad parse
							}
						}
						str = tis.readLine();
					}
					tis.close();
					
					Preconditions.checkNotNull(firstDate);
					if (lastDate == null) {
						System.out.println(file.getName()+" didn't finish, skipping");
						continue;
					}
//					Preconditions.checkState(lastDate.getTime() > firstDate.getTime(),
//							"Dates are messed up. First: "+df.format(firstDate)+", Last: "+df.format(lastDate));
					
					long diffMillis = lastDate.getTime() - firstDate.getTime();
					double diffSecs = diffMillis / 1000d;
					double diffMins = diffSecs / 60d;
					double diffHours = diffMins / 60d;
					while (diffHours < 0)
						diffHours += 24d;
					
					func.set((double)threads, diffHours);
				}
			}
		}
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(func);
		List<PlotCurveCharacterstics> chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
		PlotSpec spec = new PlotSpec(funcs, chars, "Threaded Completion Time", "# Threads", "Hours");
		
		GraphWindow gw = new GraphWindow(spec);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}

}