package scratch.aftershockStatisticsETAS;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.GregorianCalendar;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class ETAS_USGS_AftershockForecast {
	
	private static final boolean D = false; //debug
	private static final double[] min_mags_default = { 3d, 4d, 5d, 6d, 7d };
	private static final double fractile_lower = 0.025;
	private static final double fractile_median = 0.5;
	private static final double fractile_upper = 0.975;
	
	private enum Duration {
		ONE_DAY("1 Day", GregorianCalendar.DAY_OF_MONTH, 1),
		ONE_WEEK("1 Week", GregorianCalendar.WEEK_OF_MONTH, 1),
		ONE_MONTH("1 Month", GregorianCalendar.MONTH, 1),
		ONE_YEAR("1 Year", GregorianCalendar.YEAR, 1);
		
		private final String label;
		private final int calendarField;
		private final int calendarAmount;
		
		private Duration(String label, int calendarField, int calendarAmount) {
			this.label = label;
			this.calendarField = calendarField;
			this.calendarAmount = calendarAmount;
		}
		
		public GregorianCalendar getEndDate(GregorianCalendar startDate) {
			GregorianCalendar endDate = (GregorianCalendar) startDate.clone();
			endDate.add(calendarField, calendarAmount);
			return endDate;
		}
	}
	
	private ETAS_AftershockModel model;
	
	private GregorianCalendar eventDate;
	private GregorianCalendar startDate;
	private GregorianCalendar[] endDates;
	
	private Duration[] durations;
	private double[] minMags;
	private Table<Duration, Double, Double> numEventsLower;
	private Table<Duration, Double, Double> numEventsUpper;
	private Table<Duration, Double, Double> numEventsMedian;
	
	private Table<Duration, Double, Double> probs;
	
	public ETAS_USGS_AftershockForecast(ETAS_AftershockModel model,
			GregorianCalendar eventDate, GregorianCalendar startDate, double forecastEndTime) {
		this(model, min_mags_default, eventDate, startDate, forecastEndTime);
	}
	
	public ETAS_USGS_AftershockForecast(ETAS_AftershockModel model, double[] minMags,
			GregorianCalendar eventDate, GregorianCalendar startDate, double forecastEndTime) {
		compute(model, minMags, eventDate, startDate, forecastEndTime);
	}
	
	private static final DateFormat df = new SimpleDateFormat();
	
	private void compute(ETAS_AftershockModel model, double[] minMags,
			GregorianCalendar eventDate, GregorianCalendar startDate, double forecastEndTime) {
		Preconditions.checkArgument(minMags.length > 0);
		
		this.model = model;
		this.minMags = minMags;
		this.eventDate = eventDate;
		this.startDate = startDate;
		
		numEventsMedian = HashBasedTable.create();
		numEventsLower = HashBasedTable.create();
		numEventsUpper = HashBasedTable.create();
		probs = HashBasedTable.create();
		
		durations = Duration.values();
		
		GregorianCalendar forecastEndDate = new GregorianCalendar();
		forecastEndDate.setTimeInMillis(
				(long) (startDate.getTimeInMillis() + forecastEndTime*ETAS_StatsCalc.MILLISEC_PER_DAY));
		
		
		endDates = new GregorianCalendar[durations.length];
		
		double[] calcFractiles = {fractile_lower, fractile_median, fractile_upper};
		
		if(D) System.out.println("Start date: "+df.format(startDate.getTime()));
		
		for (int i=0; i<durations.length; i++) {
			Duration duration = durations[i];
						
			GregorianCalendar endDate = duration.getEndDate(startDate);
			
			if (!endDate.after(forecastEndDate)){ //check to make sure forecast extends that far
				if(D) System.out.println(duration.label+" end date: "+df.format(endDate.getTime()));

				double tMinDays = getDateDelta(eventDate, startDate);
				Preconditions.checkState(tMinDays >= 0d, "tMinDays must be positive: %s", tMinDays);
				double tMaxDays = getDateDelta(eventDate, endDate);
				Preconditions.checkState(tMaxDays > tMinDays,
						"tMaxDays must be greter than tMinDays: %s <= %s", tMaxDays, tMinDays);

				endDates[i] = endDate;


				for (int m=0; m<minMags.length; m++) {
					double minMag = minMags[m];

					double[] fractiles = model.getCumNumFractileWithAleatory(calcFractiles, minMag, tMinDays, tMaxDays);

					numEventsLower.put(duration, minMag, fractiles[0]);
					numEventsMedian.put(duration, minMag, fractiles[1]);
					numEventsUpper.put(duration, minMag, fractiles[2]);


					//				double rate = model.getModalNumEvents(minMag, tMinDays, tMaxDays);


					ArbDiscrEmpiricalDistFunc subDist = model.computeNum_DistributionFunc(tMinDays, tMaxDays, minMag);
					double prob, probMult = 1;
					if(Math.abs(subDist.getX(0) - 0) >  1e-6)
						prob = 1;
					else
						prob = 1 - subDist.getY(0);

					if(minMag < model.magComplete)
						//scale up the probability to account for events below the simulation min magnitude
						probMult = Math.pow(10, -model.b*(minMag - model.magComplete));
					prob = 1 - Math.pow(1-prob, probMult); //assumes probabilities are Poissonian

					//				System.out.println(minMag +" "+ model.magComplete +" "+ prob +" "+ probMult);
					probs.put(duration, minMag, prob);
				}
			} else {
				System.out.println("Skipping " + durations[i] + " forecast.");
				if(D) System.out.println("Forecast end date is " + df.format(endDate.getTime())
						+" but computation end date is "+ df.format(forecastEndDate.getTime()));
			}
		}
	}
	
	private static double getDateDelta(GregorianCalendar start, GregorianCalendar end) {
		long diff = end.getTimeInMillis() - start.getTimeInMillis();
		return (double)diff/(double)ProbabilityModelsCalc.MILLISEC_PER_DAY;
	}
	
	private static String[] headers = {"Forecast Interval", "Magnitude Range",
			"Median Number", "95% confidence range", "Chance of at least one"};
	
	public TableModel getTableModel() {
		final int numEach = minMags.length+1;
		final int rows = probs.rowKeySet().size()*numEach;
//		final int cols = 4;
		final int cols = 5;
		return new AbstractTableModel() {

			@Override
			public int getRowCount() {
				return rows;
			}

			@Override
			public int getColumnCount() {
				return cols;
			}

			@Override
			public Object getValueAt(int rowIndex, int columnIndex) {
				if (rowIndex == 0)
					return headers[columnIndex];
				rowIndex -= 1;
				int m = rowIndex % numEach;
				int d = rowIndex / numEach;
				
				if (columnIndex == 0) {
					if (m == 0)
						return durations[d].label;
					else if (m == 1)
						return df.format(startDate.getTime());
					else if (m == 2)
						return "through";
					else if (m == 3)
						return df.format(endDates[d].getTime());
					else
						return "";
				} else {
					if (m >= minMags.length)
						return "";
					double mag = minMags[m];
					if (columnIndex == 1) {
						return "M ≥ "+(float)mag;
					} else if (columnIndex == 2) {
						int median = (int)(numEventsMedian.get(durations[d], mag)+0.5);
						return median;
					} else if (columnIndex == 3) {
						int lower = (int)(numEventsLower.get(durations[d], mag)+0.5);
						int upper = (int)(numEventsUpper.get(durations[d], mag)+0.5);
						if (upper == 0)
							return "0";
						return lower+" to "+upper;
					} else if (columnIndex == 4) {
						int prob = (int)(100d*probs.get(durations[d], mag) + 0.5);
						if (prob == 0)
							return "<0.1%";
						else if (prob == 100)
							return ">99 %";
						return prob+" %";
					}
				}
				
				return "";
			}
			
		};
	}

}
