package scratch.aftershockStatistics;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.GregorianCalendar;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class USGS_AftershockForecast {
	
	private static final double[] min_mags_default = { 3d, 5d, 6d, 7d };
	private static final double fractile_lower = 0.025;
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
	
	private RJ_AftershockModel model;
	
	private GregorianCalendar eventDate;
	private GregorianCalendar startDate;
	private GregorianCalendar[] endDates;
	
	private Duration[] durations;
	private double[] minMags;
	private Table<Duration, Double, Double> numEventsLower;
	private Table<Duration, Double, Double> numEventsUpper;
	private Table<Duration, Double, Double> probs;
	
	public USGS_AftershockForecast(RJ_AftershockModel model,
			GregorianCalendar eventDate, GregorianCalendar startDate) {
		this(model, min_mags_default, eventDate, startDate);
	}
	
	public USGS_AftershockForecast(RJ_AftershockModel model, double[] minMags,
			GregorianCalendar eventDate, GregorianCalendar startDate) {
		compute(model, minMags, eventDate, startDate);
	}
	
	private static final DateFormat df = new SimpleDateFormat();
	
	private void compute(RJ_AftershockModel model, double[] minMags,
			GregorianCalendar eventDate, GregorianCalendar startDate) {
		Preconditions.checkArgument(minMags.length > 0);
		
		this.model = model;
		this.minMags = minMags;
		this.eventDate = eventDate;
		this.startDate = startDate;
		
		numEventsLower = HashBasedTable.create();
		numEventsUpper = HashBasedTable.create();
		probs = HashBasedTable.create();
		
		durations = Duration.values();
		endDates = new GregorianCalendar[durations.length];
		
		double[] calcFractiles = {fractile_lower, fractile_upper};
		
		System.out.println("Start date: "+df.format(startDate.getTime()));
		for (int i=0; i<durations.length; i++) {
			Duration duration = durations[i];
			GregorianCalendar endDate = duration.getEndDate(startDate);
			System.out.println(duration.label+" end date: "+df.format(endDate.getTime()));
			
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
				numEventsUpper.put(duration, minMag, fractiles[1]);
//				double rate = model.getModalNumEvents(minMag, tMinDays, tMaxDays);
				double expectedVal = model.getModalNumEvents(minMag, tMinDays, tMaxDays);
				double poissonProb = 1 - Math.exp(-expectedVal);
				probs.put(duration, minMag, poissonProb);
			}
		}
	}
	
	static double getDateDelta(GregorianCalendar start, GregorianCalendar end) {
		long diff = end.getTimeInMillis() - start.getTimeInMillis();
		return (double)diff/(double)ProbabilityModelsCalc.MILLISEC_PER_DAY;
	}
	
	private static String[] headers = {"Time Window For Analysis", "Magnitude Range",
			"Most Likely Number Of Aftershocks (95% confidence)", "Probability of one or more aftershocks"};
	
	public TableModel getTableModel() {
		final int numEach = minMags.length+1;
		final int rows = probs.rowKeySet().size()*numEach;
		final int cols = 4;
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
						int lower = (int)(numEventsLower.get(durations[d], mag)+0.5);
						int upper = (int)(numEventsUpper.get(durations[d], mag)+0.5);
						if (upper == 0)
							return "*";
						return lower+" to "+upper;
					} else if (columnIndex == 3) {
						int prob = (int)(100d*probs.get(durations[d], mag) + 0.5);
						if (prob == 0)
							return "*";
						else if (prob == 100)
							return ">99 %";
						return prob+" %";
					}
				}
				
				return "";
			}
			
		};
	}
	
	@SuppressWarnings("unchecked")
	public JSONObject buildJSON() {
		JSONObject json = new JSONObject();
		
		json.put("creationTime", System.currentTimeMillis()); // TODO is this what we want here?
		long maxEndDate = 0l;
		for (GregorianCalendar endDate : endDates)
			if (endDate.getTimeInMillis() > maxEndDate)
				maxEndDate = endDate.getTimeInMillis();
		json.put("expireTime", maxEndDate);

		// MODEL
		JSONObject modelJSON = new JSONObject();
		String name = "Reasenberg-Jones (1989, 1994) aftershock model";
		if (model instanceof RJ_AftershockModel_Bayesian)
			name += " (Bayesian Combination)";
		else if (model instanceof RJ_AftershockModel_Generic)
			name += " (Generic)";
		else if (model instanceof RJ_AftershockModel_SequenceSpecific)
			name += " (Sequence Specific)";
		modelJSON.put("name", name);
		modelJSON.put("reference", "#url");
		JSONObject modelParams = new JSONObject();
		// return AftershockStatsCalc.getExpectedNumEvents(getMaxLikelihood_a(), b, magMain, magMin, getMaxLikelihood_p(), getMaxLikelihood_c(), tMinDays, tMaxDays);
		modelParams.put("a", model.getMaxLikelihood_a());
		modelParams.put("b", model.get_b());
		modelParams.put("magMain", model.getMainShockMag());
		modelParams.put("p", model.getMaxLikelihood_p());
		modelParams.put("c", model.getMaxLikelihood_c());
		modelJSON.put("parameters", modelParams);
		json.put("model", modelJSON);
		
		// FORECAST
		JSONArray forecastsJSON = new JSONArray();
		for (int i=0; i<durations.length; i++) {
			JSONObject forecastJSON = new JSONObject();
			
			forecastJSON.put("timeStart", startDate.getTimeInMillis());
			forecastJSON.put("timeEnd", endDates[i].getTimeInMillis());
			forecastJSON.put("label", durations[i].label);
			
			JSONArray magBins = new JSONArray();
			for (int m=0; m<minMags.length; m++) {
				JSONObject magBin = new JSONObject();
				magBin.put("magnitude", minMags[m]);
				magBin.put("p95minimum", numEventsLower.get(durations[i], minMags[m]));
				magBin.put("p95maximum", numEventsUpper.get(durations[i], minMags[m]));
				magBin.put("probability", probs.get(durations[i], minMags[m]));
				magBins.add(magBin);
			}
			
			forecastJSON.put("bins", magBins);
			
			forecastsJSON.add(forecastJSON);
		}
		json.put("forecast", forecastsJSON);
		
		return json;
	}

}
