package scratch.aftershockStatistics;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.TimeZone;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class USGS_AftershockForecast {
	
	private static final double[] min_mags_default = { 3d, 5d, 6d, 7d };
	private static final double fractile_lower = 0.025;
	private static final double fractile_upper = 0.975;
	
	public enum Duration {
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
		
		@Override
		public String toString() {
			return label;
		}
	}
	
	public enum Template {
		MAINSOCK("Mainshock"),
		EQ_OF_INTEREST("Earthquake of Interest"),
		SWARM("Swarm");
		
		private String title;

		private Template(String title) {
			this.title = title;
		}
		
		@Override
		public String toString() {
			return title;
		}
	}
	
	private RJ_AftershockModel model;
	
	private int[] aftershockCounts;
	
	private GregorianCalendar eventDate;
	private GregorianCalendar startDate;
	private GregorianCalendar[] endDates;
	
	private Duration[] durations;
	private Duration advisoryDuration; // duration for the advisory paragraph
	private double[] minMags;
	private Table<Duration, Double, Double> numEventsLower;
	private Table<Duration, Double, Double> numEventsUpper;
	private Table<Duration, Double, Double> probs;
	
	private boolean includeProbAboveMainshock;
	
	// custom text which can be added
	private String injectableText = null;
	
	private Template template = Template.MAINSOCK;
	
	public USGS_AftershockForecast(RJ_AftershockModel model, List<ObsEqkRupture> aftershocks,
			GregorianCalendar eventDate, GregorianCalendar startDate) {
		this(model, aftershocks, min_mags_default, eventDate, startDate, true);
	}
	
	public USGS_AftershockForecast(RJ_AftershockModel model, List<ObsEqkRupture> aftershocks, double[] minMags,
			GregorianCalendar eventDate, GregorianCalendar startDate, boolean includeProbAboveMainshock) {
		compute(model, aftershocks, minMags, eventDate, startDate, includeProbAboveMainshock);
	}
	
	private static final DateFormat df = new SimpleDateFormat();
	private static final TimeZone utc = TimeZone.getTimeZone("UTC");
	
	private void compute(RJ_AftershockModel model, List<ObsEqkRupture> aftershocks, double[] minMags,
			GregorianCalendar eventDate, GregorianCalendar startDate, boolean includeProbAboveMainshock) {
		Preconditions.checkArgument(minMags.length > 0);
		
		this.model = model;
		this.minMags = minMags;
		this.eventDate = eventDate;
		this.startDate = startDate;
		this.includeProbAboveMainshock = includeProbAboveMainshock;
		
		// calcualte number of observations for each bin
		aftershockCounts = new int[minMags.length];
		for (int m=0; m<minMags.length; m++) {
			for (ObsEqkRupture eq : aftershocks)
				if (eq.getMag() >= minMags[m])
					aftershockCounts[m]++;
		}
		
		numEventsLower = HashBasedTable.create();
		numEventsUpper = HashBasedTable.create();
		probs = HashBasedTable.create();
		
		durations = Duration.values();
		advisoryDuration = Duration.ONE_WEEK;
		endDates = new GregorianCalendar[durations.length];
		
		double[] calcFractiles = {fractile_lower, fractile_upper};
		
		double[] calcMags = minMags;
		if (includeProbAboveMainshock) {
			// also calculate for mainshock mag
			calcMags = Arrays.copyOf(minMags, minMags.length+1);
			calcMags[calcMags.length-1] = model.getMainShockMag();
		}
		
		df.setTimeZone(utc);
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
			
			for (int m=0; m<calcMags.length; m++) {
				double minMag = calcMags[m];
				
				double[] fractiles = model.getCumNumFractileWithAleatory(calcFractiles, minMag, tMinDays, tMaxDays);
				
				numEventsLower.put(duration, minMag, fractiles[0]);
				numEventsUpper.put(duration, minMag, fractiles[1]);
//				double rate = model.getModalNumEvents(minMag, tMinDays, tMaxDays);

//				double expectedVal = model.getModalNumEvents(minMag, tMinDays, tMaxDays);
//				double poissonProb = 1 - Math.exp(-expectedVal);
				double poissonProb = model.getProbOneOrMoreEvents(minMag, tMinDays, tMaxDays);

				if (poissonProb < 1.0e-12) {
					poissonProb = 0.0;	// fewer than 4 significant digits available
				} else {
					poissonProb = Double.parseDouble (String.format ("%.3e", poissonProb));	// limit to 4 significant digits
				}

				probs.put(duration, minMag, poissonProb);
			}
		}
	}
	
	static double getDateDelta(GregorianCalendar start, GregorianCalendar end) {
		long diff = end.getTimeInMillis() - start.getTimeInMillis();
		return (double)diff/(double)ProbabilityModelsCalc.MILLISEC_PER_DAY;
	}
	
	public void setAdvisoryDuration(Duration advisoryDuration) {
		this.advisoryDuration = advisoryDuration;
	}
	
	public void setInjectableText(String injectableText) {
		this.injectableText = injectableText;
	}
	
	public boolean isIncludeProbAboveMainshock() {
		return includeProbAboveMainshock;
	}

	public void setIncludeProbAboveMainshock(boolean includeProbAboveMainshock) {
		this.includeProbAboveMainshock = includeProbAboveMainshock;
	}

	public Template getTemplate() {
		return template;
	}

	public void setTemplate(Template template) {
		this.template = template;
	}

	public Duration getAdvisoryDuration() {
		return advisoryDuration;
	}

	public String getInjectableText() {
		return injectableText;
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
	public JSONObject buildJSON () {
		return buildJSON (System.currentTimeMillis());
	}
	
	@SuppressWarnings("unchecked")
	public JSONObject buildJSON (long creation_time) {
		JSONObject json = new JSONObject();
		
		json.put("creationTime", creation_time);
		long maxEndDate = 0l;
		for (GregorianCalendar endDate : endDates)
			if (endDate.getTimeInMillis() > maxEndDate)
				maxEndDate = endDate.getTimeInMillis();
		json.put("expireTime", maxEndDate);
		
		json.put("advisoryTimeFrame", advisoryDuration.label);
		
		json.put("template", template.title);
		
		String injectableText = this.injectableText;
		if (injectableText == null)
			injectableText = "";
		json.put("injectableText", injectableText);
		
		// OBSERVATIONS
		JSONObject obsJSON = new JSONObject();
		JSONArray obsMagBins = new JSONArray();
		for (int m=0; m<minMags.length; m++) {
			JSONObject magBin = new JSONObject();
			magBin.put("magnitude", minMags[m]);
			magBin.put("count", aftershockCounts[m]);
			obsMagBins.add(magBin);
		}
		obsJSON.put("bins", obsMagBins);
		json.put("observations", obsMagBins);

		// MODEL
		JSONObject modelJSON = new JSONObject();
		modelJSON.put("name", model.getModelName());
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
				magBin.put("p95minimum", Math.round(numEventsLower.get(durations[i], minMags[m])));
				magBin.put("p95maximum", Math.round(numEventsUpper.get(durations[i], minMags[m])));
				magBin.put("probability", probs.get(durations[i], minMags[m]));
				magBins.add(magBin);
			}
			
			forecastJSON.put("bins", magBins);
			
			if (includeProbAboveMainshock) {
				JSONObject magBin = new JSONObject();
				double mainMag = model.getMainShockMag();
				magBin.put("magnitude", mainMag);
				magBin.put("probability", probs.get(durations[i], mainMag));
				forecastJSON.put("aboveMainshockMag", magBin);
			}
			
			forecastsJSON.add(forecastJSON);
		}
		json.put("forecast", forecastsJSON);
		
		return json;
	}

}
