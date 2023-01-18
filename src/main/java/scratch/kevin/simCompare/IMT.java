package scratch.kevin.simCompare;

import org.opensha.sha.imr.IntensityMeasureRelationship;
import org.opensha.sha.imr.param.IntensityMeasureParams.DurationTimeInterval;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SignificantDurationParam;

public enum IMT {
	PGV(PGV_Param.NAME, PGV_Param.NAME, PGV_Param.NAME, "pgv", "cm/s", -1d),
	PGA(PGA_Param.NAME, PGA_Param.NAME, PGA_Param.NAME, "pga", "g", 0d),
	SA0P1(SA_Param.NAME, "0.1s SA", "0.1s", "0.1s", "g", 0.1d),
	SA0P2(SA_Param.NAME, "0.2s SA", "0.2s", "0.2s", "g", 0.2d),
	SA0P5(SA_Param.NAME, "0.5s SA", "0.5s", "0.5s", "g", 0.5d),
	SA1P0(SA_Param.NAME, "1s SA", "1s", "1s", "g", 1d),
	SA2P0(SA_Param.NAME, "2s SA", "2s", "2s", "g", 2d),
	SA3P0(SA_Param.NAME, "3s SA", "3s", "3s", "g", 3d),
	SA4P0(SA_Param.NAME, "4s SA", "4s", "4s", "g", 4d),
	SA5P0(SA_Param.NAME, "5s SA", "5s", "5s", "g", 5d),
	SA6P0(SA_Param.NAME, "6s SA", "6s", "6s", "g", 6d),
	SA7P5(SA_Param.NAME, "7.5s SA", "7.5s", "7.5s", "g", 7.5),
	SA10P0(SA_Param.NAME, "10s SA", "10s", "10s", "g", 10d),
	DUR_5_75(SignificantDurationParam.NAME, "Significant Duration (5-75%)",
			"Dur5-75", "dur_5_75", "s", DurationTimeInterval.INTERVAL_5_75),
	DUR_5_95(SignificantDurationParam.NAME, "Significant Duration (5-95%)",
			"Dur5-95", "dur_5_95", "s", DurationTimeInterval.INTERVAL_5_95),
	DUR_20_80(SignificantDurationParam.NAME, "Significant Duration (20-80%)",
			"Dur20-80", "dur_20_80", "s", DurationTimeInterval.INTERVAL_20_80);
	
	private String paramName;
	private String displayName;
	private String shortName;
	private String prefix;
	private String units;
	private double period;
	private DurationTimeInterval interval;

	private IMT(String paramName, String displayName, String shortName, String prefix, String units) {
		this(paramName, displayName, shortName, prefix, units, Double.NaN);
	}
	
	private IMT(String paramName, String displayName, String shortName, String prefix, String units,
			DurationTimeInterval interval) {
		this.paramName = paramName;
		this.displayName = displayName;
		this.shortName = shortName;
		this.prefix = prefix;
		this.units = units;
		this.interval = interval;
	}
	
	private IMT(String paramName, String displayName, String shortName, String prefix, String units,
			double period) {
		this.paramName = paramName;
		this.displayName = displayName;
		this.shortName = shortName;
		this.prefix = prefix;
		this.units = units;
		this.period = period;
	}
	
	public String getParamName() {
		return paramName;
	}
	
	public String getDisplayName() {
		return displayName;
	}
	
	public String getShortName() {
		return shortName;
	}
	
	public String getPrefix() {
		return prefix;
	}
	
	public String getUnits() {
		return units;
	}
	
	public double getPeriod() {
		return period;
	}
	
	public DurationTimeInterval getDurationInterval() {
		return interval;
	}
	
	public void setIMT(IntensityMeasureRelationship imr) {
		imr.setIntensityMeasure(paramName);
		if (paramName.equals(SA_Param.NAME))
			SA_Param.setPeriodInSA_Param(imr.getIntensityMeasure(), period);
		else if (paramName.equals(SignificantDurationParam.NAME))
			SignificantDurationParam.setTimeInterval(imr.getIntensityMeasure(), interval);
	}
	
	public static IMT forPeriod(double period) {
		for (IMT imt : values())
			if (period == imt.period)
				return imt;
		throw new IllegalStateException("Period "+(float)period+" not supported");
	}
	
	public static IMT[] forPeriods(double... periods) {
		IMT[] ret = new IMT[periods.length];
		for (int p=0; p<periods.length; p++)
			ret[p] = forPeriod(periods[p]);
		return ret;
	}
	
	public static IMT forString(String string) {
		for (IMT imt : values())
			if (string.equals(imt.name()) || string.equals(imt.getDisplayName()) || string.equals(imt.getShortName()))
				return imt;
		// must be a period
		return forPeriod(Double.parseDouble(string));
	}
}