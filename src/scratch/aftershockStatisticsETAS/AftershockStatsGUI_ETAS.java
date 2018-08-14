package scratch.aftershockStatisticsETAS;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Dialog.ModalityType;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.EnumSet;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.TimeZone;
import java.util.concurrent.TimeUnit;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.data.Range;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DatasetBinner;
import org.opensha.commons.data.siteData.impl.TectonicRegime;
import org.opensha.commons.data.siteData.impl.WaldAllenGlobalVs30;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GeoTools;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.ConsoleWindow;
import org.opensha.commons.gui.plot.GraphPanel;
import org.opensha.commons.gui.plot.GraphWidget;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.param.AbstractParameter;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.param.constraint.ParameterConstraint;
import org.opensha.commons.param.constraint.impl.DoubleConstraint;
import org.opensha.commons.param.constraint.impl.RangeConstraint;
import org.opensha.commons.param.editor.AbstractParameterEditor;
import org.opensha.commons.param.editor.impl.NumericTextField;
import org.opensha.commons.param.editor.impl.ParameterListEditor;
import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.event.ParameterChangeListener;
import org.opensha.commons.param.impl.BooleanParameter;
import org.opensha.commons.param.impl.ButtonParameter;
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.commons.param.impl.EnumParameter;
import org.opensha.commons.param.impl.IntegerParameter;
import org.opensha.commons.param.impl.ParameterListParameter;
import org.opensha.commons.param.impl.RangeParameter;
import org.opensha.commons.param.impl.StringParameter;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.Blender;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.nshmp2.util.FocalMech;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupListCalc;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.CalcProgressBar;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.calc.Wald_MMI_Calc;
import org.opensha.sha.imr.param.IntensityMeasureParams.MMI_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

import scratch.aftershockStatistics.TectonicRegimeTable;
import scratch.aftershockStatistics.util.SphLatLon;
import scratch.aftershockStatistics.util.SphRegion;
import scratch.aftershockStatistics.util.SphRegionCircle;
import wContour.Global.PointD;
import wContour.Global.PolyLine;

import java.awt.Font;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;


public class AftershockStatsGUI_ETAS extends JFrame implements ParameterChangeListener {
		
	private GregorianCalendar expirationDate = new GregorianCalendar(2018, 12, 31);
	
	public AftershockStatsGUI_ETAS(String... args) {
		checkArguments(args);
		if (D) System.out.println("verbose = " + verbose + ", debug = " + D);
		if (devMode) System.out.println("Warning: Running in development mode."
				+ "This mode is untested and will probably give a bad forecast if not crash outright.");
		
		if (!devMode) {
			//check for expired software
			SimpleDateFormat formatter=new SimpleDateFormat("d MMM yyyy");
			formatter.setTimeZone(utc); //utc=TimeZone.getTimeZone("UTC"));
			double elapsedDays = (double) (System.currentTimeMillis() - expirationDate.getTimeInMillis())/ETAS_StatsCalc.MILLISEC_PER_DAY;
			
			if (elapsedDays > 0) {
				String message = "The Beta version expired on " + formatter.format(expirationDate.getTime()) + ".\n Go get the latest version from www.opensha.org/apps.";
				System.out.println(message);
				JOptionPane.showMessageDialog(null, message, "Beta version expired", JOptionPane.ERROR_MESSAGE);
				return;
			}
		}
		
		createAndShowGUI();
	}
    
	private boolean D; //debug
	private boolean devMode;
	private boolean verbose;
	private volatile boolean changeListenerEnabled = true;
	private boolean tipsOn = true;
	private File workingDir;

	//this is needed to prevent long processing times/overloaded memory. If more than MAX_EARTHQUAKE_NUMBER aftershocks are retrieved from ComCat,
	//the magnitude of completeness is increased to the size of the nth largest aftershock. Must be > 0  or you'll get a nullPoitnerException down the line.
	private final static int MAX_EARTHQUAKE_NUMBER = 1000; 
	
	
	/*
	 * Data parameters
	 */
	private static final long serialVersionUID = 1L;
	private StringParameter eventIDParam;
	
	private ParameterList timeWindow;
	private ParameterListParameter timeWindowEditParam;
	private ButtonParameter quickForecastButton;
	private DoubleParameter dataStartTimeParam;
	private DoubleParameter dataEndTimeParam;
	private BooleanParameter nowBoolean;
	
	private JFileChooser workingDirChooser;
	
	private enum ForecastDuration {
		
		YEAR("year",365d),
		MONTH("month",30d),
		WEEK("week", 7d),
		DAY("day",1d);
		
		private double duration;
		private String name;
		
		private ForecastDuration(String name, double duration) {
			this.duration = duration;
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
	}
	private BooleanParameter plotAllDurationsParam;
	
	private enum RegionType {
		CIRCULAR("Circular"),
		CIRCULAR_WC94("Automatic"),
		RECTANGULAR("Rectangular");
		
		private String name;
		
		private RegionType(String name) {
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
		
		public boolean isCircular() {
			return this == CIRCULAR_WC94 || this == CIRCULAR;
		}
	}
	
	private enum RegionCenterType {
		EPICENTER("Epicenter"),
		CENTROID("Aftershock Centroid"),
		SPECIFIED("Custom Location");
		
		private String name;
		
		private RegionCenterType(String name) {
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
	}
	
	
	
	private EnumParameter<RegionType> regionTypeParam;
	private final static String DECIMAL_DEGREES = "Decimal Degrees";
	
	private DoubleParameter regionCenterLatParam;
	private DoubleParameter regionCenterLonParam;
	
	private DoubleParameter radiusParam;
	private DoubleParameter minLatParam;
	private DoubleParameter maxLatParam;
	private DoubleParameter minLonParam;
	private DoubleParameter maxLonParam;
	private DoubleParameter minDepthParam;
	private DoubleParameter maxDepthParam;
	private EnumParameter<RegionCenterType> regionCenterTypeParam;
	
	private ParameterList constraintList;
	private ParameterListParameter constraintEditParam;
	
	private ParameterList regionList;
	private ParameterListParameter regionEditParam;
	private EnumParameter<TectonicRegime> tectonicRegimeParam;
		
	private ButtonParameter fetchButton;
	
	private JFileChooser loadCatalogChooser;
	private ButtonParameter loadCatalogButton;
	
	private JFileChooser saveCatalogChooser;
	private ButtonParameter saveCatalogButton;
	
	/*
	 * B-value fit parameters
	 */
	private BooleanParameter autoMcParam;
	private DoubleParameter mcParam;
	private DoubleParameter magPrecisionParam;
	private ButtonParameter computeBButton;
	
	
	/*
	 * Aftershock model parameters
	 */
	
	private DoubleParameter magRefParam;
	private BooleanParameter fitMSProductivityParam;
	
	private RangeParameter amsValRangeParam;
	private IntegerParameter amsValNumParam;
	private RangeParameter aValRangeParam;
	private IntegerParameter aValNumParam;
	private RangeParameter pValRangeParam;
	private IntegerParameter pValNumParam;
	private RangeParameter cValRangeParam;
	private IntegerParameter cValNumParam;
	
	private BooleanParameter timeDepMcParam;
	private DoubleParameter rmaxParam;

	private ButtonParameter computeAftershockParamsButton;
	
	private DoubleParameter amsValParam;
	private DoubleParameter aValParam;
	private DoubleParameter pValParam;
	private DoubleParameter cValParam;
	private DoubleParameter bParam;
	private DoubleParameter alphaParam;
	
	private DoubleParameter mapLevelParam;
	private DoubleParameter mapScaleParam;	
	private DoubleParameter mapGridSpacingParam;
	private DoubleParameter mapPOEParam;
	
	private DoubleParameter forecastStartTimeParam;
	private DoubleParameter forecastEndTimeParam;
	private DoubleParameter forecastDurationParam;
	
	private ButtonParameter computeAftershockForecastButton;
	private BooleanParameter plotSpecificOnlyParam;
	private ButtonParameter generateMapButton;
	
	private BooleanParameter vs30Param;
	
	private enum IntensityType {
		MMI("MMI (Modified Mercali Intensity)", "MMI"),
		PGA("Peak acceleration (%g)", "PGA"),
		PGV("Peak velocity (cm/s)", "PGV"),
		PSA("Peak Spectral Acceleration 1Hz (%g)", "PSA"),
		NONE("Aftershock Rate Only", "NONE");
		
		private String name;
		private String abbreviation;
		
		private IntensityType(String name, String abbreviation) {
			this.name = name;
			this.abbreviation = abbreviation;
		}
		
		@Override
		public String toString() {
			return name;
		}
		
		public String getAbbreviation() {
			return this.abbreviation;
		}
	}
	private EnumParameter<IntensityType> intensityTypeParam;	
	
	private enum MapType {
		PROBABILITIES("Probabilities"),
		LEVEL("Levels");
		
		private String name;
		
		private MapType(String name) {
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
	}
	private EnumParameter<MapType> mapTypeParam;	
	
	private enum FitSourceType {
		SHAKEMAP("Shakemap finite source"),
		AFTERSHOCKS("Early aftershocks"),
		POINT_SOURCE("Point source"),
		CUSTOM("Load from file");
		
		private String name;
		
		private FitSourceType(String name) {
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
	}
	private EnumParameter<FitSourceType> fitSourceTypeParam;	
	
	private ButtonParameter publishAdvisoryButton;
	
	private JTabbedPane tabbedPane;
	private JScrollPane consoleScroll;
	
	private final int console_tab_index = 0;
	private final int catalog_tab_index = 1;
	private final int epicenter_tab_index = 2;
	private final int mag_time_tab_index = 3;
	private final int mag_num_tab_index = 4;
	private final int cml_num_tab_index = 5;
	private final int pdf_tab_index = 6;
	private final int aftershock_expected_index = 7;
	private final int forecast_table_tab_index = 8;
	private final int forecast_map_tab_index = 9;
	
	private GraphWidget epicenterGraph;
	private GraphWidget magNumGraph;
	private GraphWidget magTimeGraph;
	private GraphWidget cmlNumGraph;
	private JScrollPane catalogPane;
	private JTextArea catalogText;
	private JTabbedPane pdfGraphsPane;
	
	private List<GraphWidget> aftershockExpectedNumGraph;
	private List<GraphWidget> aftershockProbabilityGraph;
	
	private JTabbedPane forecastTablePane;
	private JTabbedPane forecastMFDPane;
	private JTabbedPane forecastMapPane;
	
	private ParameterList forecastParams;
	private ParameterList dataParams;
	private ParameterList mfdParams;
	private ParameterList fitParams;
	private ParameterList mapParams;
	private ParameterList outputParams;
	private ParameterList mapPlotParams;
	private ParameterList publishAdvisoryParams;
	
	private ParameterListEditor forecastEditor;
	private IncrementalMagFreqDist aftershockMND;
	private ParameterListEditor dataEditor;
	private ParameterListEditor mfdEditor;
	private ParameterListEditor fitEditor;
	private ParameterListEditor mapEditor;
	private ParameterListEditor outputEditor;
	private ParameterListEditor mapPlotEditor;
	private ParameterListEditor publishAdvisoryEditor;
	
	private ETAS_ComcatAccessor accessor;
	private WC1994_MagLengthRelationship wcMagLen;
	
	private SphRegion region;
	private ObsEqkRupture mainshock;
	private FaultTrace faultTrace;
	private ObsEqkRupList aftershocks;
	private ObsEqkRupture largestShock;
	private ETAS_RateModel2D rateModel2D;
	private List<ContourModel> contourList;
	private String shakeMapURL;
	private DiscretizedFunc[] pgvCurves = null; //defining this globally so it doesn't need to be recomputed
	private DiscretizedFunc[] pgaCurves = null; //defining this globally so it doesn't need to be recomputed
	private DiscretizedFunc[] psaCurves = null; //defining this globally so it doesn't need to be recomputed
	
	private OgataMagFreqDist ogataMND;
		
	private GenericETAS_ParametersFetch genericFetch = null;
	private GenericETAS_Parameters genericParams = null;
	private ETAS_AftershockModel_Generic genericModel = null;
//	private ETAS_AftershockModel_SequenceSpecific seqSpecModel = null;
//	private ETAS_AftershockModel_Bayesian bayesianModel = null; // the bayesian model is now just a sequence specific model
	private ETAS_AftershockModel_SequenceSpecific bayesianModel = null;
	
	private final Color generic_color = new Color(136,119,156);
	private final Color bayesian_color = new Color(84,185,157);
	private final Color sequence_specific_color = new Color(31,163,218);
		
	
	private int height = 750;
	private int paramWidth = 220;
	private int outputWidth = 65;
	private int chartWidth = 750;
	private int chartHeight = height;
	private int sigDigits = 4;
	
	/* Initialize the GUI
	 *  
	 */
	private void createAndShowGUI() {
		/*
		 * Data parameters
		 */
		forecastParams = new ParameterList();
		dataParams = new ParameterList();
		mfdParams = new ParameterList();
		fitParams = new ParameterList();
		mapParams = new ParameterList();
		outputParams = new ParameterList();
		mapPlotParams = new ParameterList();
		publishAdvisoryParams = new ParameterList();
		
		eventIDParam = new StringParameter("USGS Event ID", "");
//		eventIDParam.setValue(null);
		eventIDParam.setInfo("the event ID can be found in the event page URL for specific earthquakes at https://earthquake.usgs.gov/earthquakes/");
		eventIDParam.addParameterChangeListener(this);
		forecastParams.addParameter(eventIDParam);
		
		nowBoolean = new BooleanParameter("Begin at current time", true);
		nowBoolean.setInfo("Use current time for forecast start");
		nowBoolean.addParameterChangeListener(this);
		forecastParams.addParameter(nowBoolean);
		
		forecastStartTimeParam = new DoubleParameter("Forecast Start Time", 0d, Double.POSITIVE_INFINITY);
		forecastStartTimeParam.setUnits("Days");
		forecastStartTimeParam.setInfo("Time at which to begin forecast, relative to main shock origin time");
		forecastStartTimeParam.addParameterChangeListener(this);
		forecastParams.addParameter(forecastStartTimeParam);

		forecastDurationParam = new DoubleParameter("Forecast Duration", 0d, 366, new Double(366d));
		forecastDurationParam.setUnits("Days");
		forecastDurationParam.setInfo("Forecast duration realative to forecast start time");
		forecastDurationParam.addParameterChangeListener(this);
		forecastParams.addParameter(forecastDurationParam);

		plotAllDurationsParam = new BooleanParameter("Plot Week/Month/Year", true);	
		plotAllDurationsParam.addParameterChangeListener(this);
		plotAllDurationsParam.setInfo("Build plots for day, week, month, and year increments");
		forecastParams.addParameter(plotAllDurationsParam);
		
		quickForecastButton = new ButtonParameter("Forecast using default settings", "Quick Forecast");
		quickForecastButton.setInfo("Run all steps using default settings");
		quickForecastButton.addParameterChangeListener(this);
		quickForecastButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		forecastParams.addParameter(quickForecastButton);
		
		timeWindow = new ParameterList(); 
		dataStartTimeParam = new DoubleParameter("Data Start Time", 0d, 366, new Double(0d));
		dataStartTimeParam.setUnits("Days");
		dataStartTimeParam.setInfo("Relative to main shock origin time");
		dataStartTimeParam.addParameterChangeListener(this);
		timeWindow.addParameter(dataStartTimeParam);
		
		dataEndTimeParam = new DoubleParameter("Data End Time", 0d, Double.POSITIVE_INFINITY);
		dataEndTimeParam.setUnits("Days");
		dataEndTimeParam.setInfo("Relative to main shock origin time");
		dataEndTimeParam.addParameterChangeListener(this);
		timeWindow.addParameter(dataEndTimeParam);
				
		forecastEndTimeParam = new DoubleParameter("Forecast End Time", 0d, Double.POSITIVE_INFINITY);
		forecastEndTimeParam.setUnits("Days");
//		forecastEndTimeParam.addParameterChangeListener(this);
//		timeWindow.addParameter(forecastEndTimeParam);
		
		timeWindowEditParam = new ParameterListParameter("Edit data time window", timeWindow);
		if (devMode) forecastParams.addParameter(timeWindowEditParam);		
		
		// these are inside region editor
		regionTypeParam = new EnumParameter<AftershockStatsGUI_ETAS.RegionType>(
				"Aftershock Zone Type", EnumSet.allOf(RegionType.class), RegionType.CIRCULAR_WC94, null);
		regionTypeParam.setInfo("Different ways to specify the aftershock zone");
		regionTypeParam.addParameterChangeListener(this);
//		dataParams.addParameter(regionTypeParam);

		regionCenterLatParam = new DoubleParameter("Center Latitude",
				new DoubleConstraint(GeoTools.LAT_MIN,GeoTools.LAT_MAX),
				DECIMAL_DEGREES, null);
		regionCenterLatParam.getConstraint().setNullAllowed(true);
		regionCenterLatParam.setInfo("Center point of aftershock zone");

		regionCenterLonParam = new DoubleParameter("Center Longitude",
				new DoubleConstraint(GeoTools.LON_MIN,360),
				DECIMAL_DEGREES, null);
		regionCenterLonParam.getConstraint().setNullAllowed(true);
		regionCenterLonParam.setInfo("Center point of aftershock zone");
		
		radiusParam = new DoubleParameter("Radius", 0d, 10000d, new Double(20));
		radiusParam.setUnits("km");
		radiusParam.setInfo("Radius for collecting aftershocks in km");
		minLatParam = new DoubleParameter("Min Lat", -90d, 90d, new Double(0));
		minLatParam.setInfo("Minimum latitude for rectangle");
		maxLatParam = new DoubleParameter("Max Lat", -90d, 90d, new Double(0));
		maxLatParam.setInfo("Maximum latitude for rectangle");
		minLonParam = new DoubleParameter("Min Lon", -180d, 360d, new Double(0));
		minLonParam.setInfo("Minimum longitude for rectangle");
		maxLonParam = new DoubleParameter("Max Lon", -180d, 360d, new Double(0));
		maxLonParam.setInfo("Maximum longitude for rectangle");
		
		minDepthParam = new DoubleParameter("Min Depth", 0d, 1000d, new Double(0));
		minDepthParam.setUnits("km");
		minDepthParam.setInfo("Minimum depth for collecting aftershocks");
		maxDepthParam = new DoubleParameter("Max Depth", 0d, 1000d, new Double(1000d));
		maxDepthParam.setUnits("km");
		maxDepthParam.setInfo("Maximum depth for collecting aftershocks");
		
		regionCenterTypeParam = new EnumParameter<AftershockStatsGUI_ETAS.RegionCenterType>(
				"Aftershock Zone Center", EnumSet.allOf(RegionCenterType.class), RegionCenterType.CENTROID, null);

//		regionCenterLocParam = new LocationParameter("Aftershock Zone Center Location");
		regionCenterTypeParam.setInfo("Different ways to select the aftershock zone center");
		regionCenterTypeParam.addParameterChangeListener(this);
		
		regionList = new ParameterList();
		
		regionEditParam = new ParameterListParameter("Edit Aftershock Zone", regionList);
		regionEditParam.setInfo("Manually edit the aftershock collection region");
		regionEditParam.addParameterChangeListener(this);
		regionEditParam.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		dataParams.addParameter(regionEditParam);
		
		fetchButton = new ButtonParameter("USGS Event Webservice", "Fetch Data");
		fetchButton.setInfo("Download aftershocks from USGS ComCat");
		fetchButton.addParameterChangeListener(this);
		fetchButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		dataParams.addParameter(fetchButton);
		
		loadCatalogButton = new ButtonParameter("External Catalog", "Load Local Catalog");
		loadCatalogButton.setInfo("Load local catalog from disk");
		loadCatalogButton.addParameterChangeListener(this);
		loadCatalogButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		dataParams.addParameter(loadCatalogButton);
		
		saveCatalogButton = new ButtonParameter("Aftershock Catalog", "Save Local Catalog");
		saveCatalogButton.setInfo("Save local catalog to disk");
		saveCatalogButton.addParameterChangeListener(this);
		saveCatalogButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		dataParams.addParameter(saveCatalogButton);
		
		
		/*
		 * Constraint params
		 */
		computeBButton = new ButtonParameter("b-value", "Compute b");
		computeBButton.setInfo("Compute sequence-specific Gutenberg-Richter b-value");
		computeBButton.addParameterChangeListener(this);
		computeBButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		mfdParams.addParameter(computeBButton);
		
		autoMcParam = new BooleanParameter("Automatically find Mc", true);
//		autoMcParam.addParameterChangeListener(this);
//		mfdParams.addParameter(autoMcParam);
		
		bParam = new DoubleParameter("b", 1d);
		bParam.setInfo("Gutenberg-Richter b-value");
		bParam.addParameterChangeListener(this);
		bParam.getEditor().getComponent().setMinimumSize(null);
		bParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		
		mcParam = new DoubleParameter("Mc", 5.0);
		mcParam.setDefaultValue(5.0);
		mcParam.setInfo("Magnitude of completeness");
		mcParam.addParameterChangeListener(this);
		mcParam.getEditor().getComponent().setMinimumSize(null);
		mcParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		mcParam.getEditor().getComponent().setAlignmentX(RIGHT_ALIGNMENT);
		
		
		magPrecisionParam = new DoubleParameter("\u0394M", 0d, 1d, new Double(0.1));
		magPrecisionParam.setInfo("Magnitude rounding applied by network");;
		magPrecisionParam.addParameterChangeListener(this);
		magPrecisionParam.getEditor().getComponent().setMinimumSize(null);
		magPrecisionParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		
		
		
		if (!devMode){
			bParam.getEditor().setEnabled(devMode);
			setEnabledStyle(bParam, false);
		
//			mcParam.getEditor().setEnabled(devMode); // leave this one editable
			magPrecisionParam.getEditor().setEnabled(devMode);
		}
		
		outputParams.addParameter(mcParam);
		outputParams.addParameter(bParam);
		
//		outputParams.addParameter(magPrecisionParam);
		
		
		
//		mfdParams.addParameter(bParam);
//		mfdParams.addParameter(mcParam);
//		mfdParams.addParameter(magPrecisionParam);
		
		alphaParam = new DoubleParameter("alpha-value", 1d);
		alphaParam.setInfo("Linked to b-value");
		alphaParam.addParameterChangeListener(this);
//		dataParams.addParameter(alphaParam);
		
		/* make this display more useful labels */
		tectonicRegimeParam = new EnumParameter<TectonicRegime>(
				"Tectonic Regime", EnumSet.allOf(TectonicRegime.class), TectonicRegime.GLOBAL_AVERAGE, null);
		tectonicRegimeParam.setInfo("Choose tectonic regime. ACR = active continental region, SCR = stable continental region, SOR = stable oceanic region, SZ = subduction zone");
		tectonicRegimeParam.addParameterChangeListener(this);
		fitParams.addParameter(tectonicRegimeParam);
		
		
		//these are inside constraint editor
		fitMSProductivityParam = new BooleanParameter("Fit MS Productivity", true);
		fitMSProductivityParam.addParameterChangeListener(this);
		
		amsValRangeParam = new RangeParameter("ams-value range", new Range(-4.0, -1));
		amsValRangeParam.setInfo("Specify mainshock productivity parameter range.");
		amsValRangeParam.addParameterChangeListener(this);
		
		amsValNumParam = new IntegerParameter("ams-value num", 1, 101, new Integer(51));
		amsValNumParam.setInfo("Set number of points in grid search");
		amsValNumParam.addParameterChangeListener(this);
		
		aValRangeParam = new RangeParameter("a-value range", new Range(-4.0, -1));
		aValRangeParam.setInfo("Specify secondary productivity parameter range.");
		aValRangeParam.addParameterChangeListener(this);
		
		aValNumParam = new IntegerParameter("a-value num", 1, 101, new Integer(31));
		aValNumParam.setInfo("Set number of points in grid search");
		aValNumParam.addParameterChangeListener(this);
		
		pValRangeParam = new RangeParameter("p-value range", new Range(0.5, 2.0));
		pValRangeParam.setInfo("Specify Omori p-parameter range.");
		pValRangeParam.addParameterChangeListener(this);
		
		pValNumParam = new IntegerParameter("p-value num", 1, 101, new Integer(31));
		pValNumParam.setInfo("Set number of points in grid search");
		pValNumParam.addParameterChangeListener(this);
		
		cValRangeParam = new RangeParameter("c-value range", new Range(1e-5, 1));
		cValRangeParam.setInfo("Specify Omori c-parameter range.");
		cValRangeParam.addParameterChangeListener(this);
		
		cValNumParam = new IntegerParameter("c-value num", 1, 101, new Integer(21));
		cValNumParam.setInfo("Set number of points in grid search");
		cValNumParam.addParameterChangeListener(this);
		
		rmaxParam = new DoubleParameter("rmax", 1d, Double.POSITIVE_INFINITY, new Double(200));
		rmaxParam.setInfo("Specify max completeness rate");
		rmaxParam.addParameterChangeListener(this);
		
		magRefParam = new DoubleParameter("Reference Magnitude", new Double(0d));
		
		timeDepMcParam = new BooleanParameter("Apply time dep. Mc", true);
		timeDepMcParam.setInfo("Apply time dependent magnitude of completeness");
		timeDepMcParam.addParameterChangeListener(this);

		constraintList = new ParameterList();
		if(devMode) constraintList.addParameter(fitMSProductivityParam);
		constraintList.addParameter(timeDepMcParam);
		constraintList.addParameter(amsValRangeParam);
		constraintList.addParameter(amsValNumParam);
		constraintList.addParameter(aValRangeParam);
		constraintList.addParameter(aValNumParam);
		constraintList.addParameter(pValRangeParam);
		constraintList.addParameter(pValNumParam);
		constraintList.addParameter(cValRangeParam);
		constraintList.addParameter(cValNumParam);
		if(devMode) constraintList.addParameter(rmaxParam);
				
		constraintEditParam = new ParameterListParameter("Edit Fit Constraints", constraintList);
		constraintEditParam.setInfo("Manually edit aftershock model constraints");
		constraintEditParam.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		constraintEditParam.addParameterChangeListener(this);
		fitParams.addParameter(constraintEditParam);
				
		computeAftershockParamsButton = new ButtonParameter("Aftershock Params", "Compute Model Fit");
		computeAftershockParamsButton.addParameterChangeListener(this);
		computeAftershockParamsButton.setInfo("Estimate aftershock model from observed aftershocks");
		computeAftershockParamsButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		fitParams.addParameter(computeAftershockParamsButton);
		
		/*
		 * Fit params
		 */
		amsValParam = new DoubleParameter("ams", new Double(0d));
		amsValParam.setInfo("Mainshock direct productivity parameter");
		amsValParam.setValue(null);
		
		amsValParam.addParameterChangeListener(this);
		amsValParam.getEditor().getComponent().setMinimumSize(null);
		amsValParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		setEnabledStyle(amsValParam, false);
		outputParams.addParameter(amsValParam);
		
		aValParam = new DoubleParameter("a", new Double(0d));
		aValParam.setInfo("Secondary productivity parameter");
		aValParam.setValue(null);
		aValParam.addParameterChangeListener(this);
		aValParam.getEditor().getComponent().setMinimumSize(null);
		aValParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		setEnabledStyle(aValParam, false);
		outputParams.addParameter(aValParam);
		
		pValParam = new DoubleParameter("p", new Double(0d));
		pValParam.setValue(null);
		pValParam.setInfo("Omori p-parameter");
		pValParam.addParameterChangeListener(this);
		pValParam.getEditor().getComponent().setMinimumSize(null);
		pValParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		setEnabledStyle(pValParam, false);
		outputParams.addParameter(pValParam);
		
		cValParam = new DoubleParameter("log-c", new Double(0d));
		cValParam.setValue(null);
		cValParam.setInfo("log10 of Omori c-parameter");
		cValParam.addParameterChangeListener(this);
		cValParam.getEditor().getComponent().setMinimumSize(null);
		cValParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		setEnabledStyle(cValParam, false);
		outputParams.addParameter(cValParam);

		computeAftershockForecastButton = new ButtonParameter("Aftershock Forecast", "Run Generic Forecast");
		computeAftershockForecastButton.setInfo("Compute aftershock forecast using typical parameters");
		computeAftershockForecastButton.addParameterChangeListener(this);
		computeAftershockForecastButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		fitParams.addParameter(computeAftershockForecastButton);
		
		plotSpecificOnlyParam = new BooleanParameter("Specific Forecast Only", false);
		plotSpecificOnlyParam.setInfo("Check to plot only the sequence-specific forecast in summary figures.");
		plotSpecificOnlyParam.addParameterChangeListener(this);
		plotSpecificOnlyParam.getEditor().setEnabled(false);
		plotSpecificOnlyParam.getEditor().refreshParamEditor();
		fitParams.addParameter(plotSpecificOnlyParam);
		
		/*
		 * map params
		 */
		intensityTypeParam = new EnumParameter<IntensityType>("Intensity Measure",
				EnumSet.allOf(IntensityType.class), IntensityType.MMI, null);
		intensityTypeParam.setInfo("Specify shaking intensity metric for map");
		intensityTypeParam.addParameterChangeListener(this);
		mapParams.addParameter(intensityTypeParam);
		
		mapTypeParam =  new EnumParameter<MapType>("Map Type",
				EnumSet.allOf(MapType.class), MapType.PROBABILITIES, null);
		mapTypeParam.setInfo("Specify map type: Probabilities of exceeding a fixed level, or Levels with a fixed probability of exceedance");
//		if (devMode) {
			mapTypeParam.addParameterChangeListener(this);
			mapParams.addParameter(mapTypeParam);
//		}
		
		fitSourceTypeParam = new EnumParameter<FitSourceType>("Fit finite source",
				EnumSet.allOf(FitSourceType.class), FitSourceType.SHAKEMAP, null);
		fitSourceTypeParam.setInfo("Options for fitting the spatial extent of the mainshock rupture");
		fitSourceTypeParam.addParameterChangeListener(this);
		mapParams.addParameter(fitSourceTypeParam);
		
		vs30Param = new BooleanParameter("Apply site corrections", true);
		vs30Param.setInfo("Adjust estimated ground motions based on local site conditions (recommended)");
		vs30Param.addParameterChangeListener(this);
		mapParams.addParameter(vs30Param);
		
		generateMapButton = new ButtonParameter("Forecast Map", "Render");
		generateMapButton.setInfo("Generate map(s)");
		generateMapButton.addParameterChangeListener(this);
		generateMapButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		mapParams.addParameter(generateMapButton);
	
		mapLevelParam = new DoubleParameter("Level", 1d, 100d, new Double(10) );
		mapLevelParam.setInfo("Shaking intensity level");
		mapLevelParam.addParameterChangeListener(this);
		mapLevelParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		mapPlotParams.addParameter(mapLevelParam);
		
		mapPOEParam = new DoubleParameter("POE (%)", 0, 99.9, new Double(10));
		mapPOEParam.setInfo("Probability of exceedence");
//		if (devMode) {
			mapPOEParam.addParameterChangeListener(this);
			mapPOEParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
			mapPlotParams.addParameter(mapPOEParam);
//		}
		
		mapGridSpacingParam = new DoubleParameter("\u0394 (km)", 0.01, 10000, new Double(10d));
		mapGridSpacingParam.setInfo("Cell size for map (km)");
		mapGridSpacingParam.addParameterChangeListener(this);
		mapGridSpacingParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		mapPlotParams.addParameter(mapGridSpacingParam);
		
		mapScaleParam  = new DoubleParameter("Scale", 0.5d, 10d, new Double(1));
		mapScaleParam.setInfo("Factor by which to extend the map region beyond the aftershock zone");
		mapScaleParam.addParameterChangeListener(this);
		mapScaleParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		mapPlotParams.addParameter(mapScaleParam);
		
		publishAdvisoryButton = new ButtonParameter("Publish Advisory", "Publish");
		publishAdvisoryButton.addParameterChangeListener(this);
		publishAdvisoryButton.setInfo("Create and save forecast summary documents");
		publishAdvisoryButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		publishAdvisoryParams.addParameter(publishAdvisoryButton);
		
		ConsoleWindow console = new ConsoleWindow(true);
		consoleScroll = console.getScrollPane();
		
		JTextArea text = console.getTextArea();
		text.setCaretPosition(0);
		text.setCaretPosition(text.getText().length());

		tabbedPane = new JTabbedPane();
		tabbedPane.setTabLayoutPolicy(JTabbedPane.WRAP_TAB_LAYOUT);

		JPanel mainPanel = new JPanel(new BorderLayout());
		JPanel displayPanel = new JPanel(new BorderLayout());
		JPanel paramsPanel = new JPanel(new BorderLayout());
		JPanel dataParamsPanel = new JPanel(new BorderLayout());
		JPanel fitParamsPanel = new JPanel(new BorderLayout());
		JPanel outputPanel = new JPanel(new BorderLayout());
		
		forecastEditor = new ParameterListEditor(forecastParams);
		dataEditor = new ParameterListEditor(dataParams);
		mfdEditor = new ParameterListEditor(mfdParams);
		fitEditor = new ParameterListEditor(fitParams);
		mapEditor = new ParameterListEditor(mapParams);
		outputEditor = new ParameterListEditor(outputParams);
		mapPlotEditor = new ParameterListEditor(mapPlotParams);
		publishAdvisoryEditor = new ParameterListEditor(publishAdvisoryParams);
		
		forecastEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/10d*(5 + 0.0))));
		dataEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/10d*(5 + 0.0))));
		
		mfdEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/10d*(1 + 0.0))));
		fitEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/10d*(4.3 + 0.0))));
		mapEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/10d*(4.7 + 0.0))));
		publishAdvisoryEditor.setPreferredSize(new Dimension(paramWidth, (int)(height/10d*(1 + 0.0))));
		
		outputEditor.setPreferredSize(new Dimension((int) (paramWidth/4d), (int) (height/11d * 7)));
		mapPlotEditor.setPreferredSize(new Dimension((int) (paramWidth/4d), (int) (height/11d * 4)));
		
		displayPanel.setPreferredSize(new Dimension(chartWidth, height));
		consoleScroll.setPreferredSize(new Dimension(chartWidth, height));
		paramsPanel.setPreferredSize(new Dimension(2*paramWidth+outputWidth+25, height));
		dataParamsPanel.setPreferredSize(new Dimension(paramWidth, height));
		fitParamsPanel.setPreferredSize(new Dimension(paramWidth, height));
		outputPanel.setPreferredSize(new Dimension(outputWidth+25, height));
		
		forecastEditor.setTitle("Forecast parameters");
		publishAdvisoryEditor.setTitle("Publish Advisory");
		dataEditor.setTitle("Fetch aftershock data");
		mfdEditor.setTitle("Fit magnitude distribution");
		fitEditor.setTitle("Calculate forecast");
		mapEditor.setTitle("Generate map");
		outputEditor.setTitle("Params");
		mapPlotEditor.setTitle("Map Opts.");
		
		dataParamsPanel.add(forecastEditor, BorderLayout.NORTH);
		dataParamsPanel.add(dataEditor, BorderLayout.CENTER);
		dataParamsPanel.add(mfdEditor, BorderLayout.SOUTH);		
		
		fitParamsPanel.add(fitEditor, BorderLayout.NORTH);
		fitParamsPanel.add(mapEditor, BorderLayout.CENTER);
		fitParamsPanel.add(publishAdvisoryEditor, BorderLayout.SOUTH);
		
		outputPanel.add(outputEditor, BorderLayout.CENTER);
		outputPanel.add(mapPlotEditor, BorderLayout.SOUTH);

		paramsPanel.add(dataParamsPanel, BorderLayout.WEST);
		paramsPanel.add(fitParamsPanel, BorderLayout.CENTER);
		paramsPanel.add(outputPanel, BorderLayout.EAST);
		
		displayPanel.add(tabbedPane, BorderLayout.CENTER);
		
		//initialize tabs (careful to do it in the same order as in the list)
		tabbedPane.addTab("Console", null, consoleScroll, "Console with status and error messages (technical)");
		
		catalogText = new JTextArea();
		catalogText.setEditable(false);
		catalogPane = new JScrollPane(catalogText);
		tabbedPane.addTab("Catalog", null, catalogPane, "Observed aftershock catalog");

		tabbedPane.addTab("Epicenters", null, epicenterGraph, "Observed aftershock epicenter map");
		
		tabbedPane.addTab("Mag/Time Plot", null, magTimeGraph,
				"Observed aftershock magnitudes vs. time");

		tabbedPane.addTab("Mag/Num Plot", null, magNumGraph,
				"Observed aftershock magnitude-number distribution");
		
		tabbedPane.addTab("Cml Number Plot", null, cmlNumGraph,
				"Observed and Forecast number of aftershocks with time");

		tabbedPane.addTab("Model PDFs", null, pdfGraphsPane,
				"Probability distributions for model parameter estimates (technical)");

		tabbedPane.addTab("Forecast MFD", null, forecastMFDPane,
				"Forecast magnitude-number and magnitude-probability distributions");

		tabbedPane.addTab("Forecast Table", null, forecastTablePane,
				"Aftershock forecast table");

		tabbedPane.addTab("Forecast Maps", null, forecastMapPane,
				"Aftershock forecast maps");

		// disable tabs beyond console tab.
		for (int i = 1; i < tabbedPane.getTabCount() ; i++){
			tabbedPane.setForegroundAt(i, new Color(128,128,128));
			tabbedPane.setEnabledAt(i, false);
		}

		mainPanel.add(paramsPanel, BorderLayout.WEST);
		mainPanel.add(displayPanel, BorderLayout.CENTER);

		tabbedPane.setSelectedIndex(0);

		setContentPane(mainPanel);
		setSize(paramWidth*2 + outputWidth + 25 + chartWidth, height);
		setVisible(true);
		
		setDefaultCloseOperation(EXIT_ON_CLOSE);
		setTitle("Aftershock Forecaster (Beta)");
		setLocationRelativeTo(null);
		
		workingDir = new File(System.getProperty("user.home"));
		accessor = new ETAS_ComcatAccessor();

		updateRegionParamList(regionTypeParam.getValue(), regionCenterTypeParam.getValue());
		updateMapPanel();
		
		setEnableParameterEditing(devMode);
		refreshTimeWindowEditor();
		setEnableParamsPostFetch(false);
		setEnableParamsPostForecast(false);
		 
		printTip(0);
	} //end createAndShowGUI()
	
	/*
	 * Sets the style of the parameter to reflect non-editable (output) style
	 */
	private void setEnabledStyle(AbstractParameter<?> param, boolean enabled) {
		if (enabled) {
			if (param.getEditor() instanceof AbstractParameterEditor) {
				Font font = param.getEditor().getComponent().getFont();
				AbstractParameterEditor<?> editor = (AbstractParameterEditor<?>) param.getEditor();
				// update title font
				editor.setTitleFont(new Font(font.getFamily(), Font.BOLD, font.getSize()));
				if (editor.getWidget() instanceof NumericTextField) {
					NumericTextField textField = (NumericTextField)editor.getWidget();
					textField.setFont(new Font(font.getFamily(), Font.PLAIN, font.getSize()));
					textField.setBackground(new Color(255,255,255));
					textField.setForeground(Color.BLACK);
			 	}
				editor.repaint();
				
			}
		} else {
			if (param.getEditor() instanceof AbstractParameterEditor) {
				Font font = param.getEditor().getComponent().getFont();
				AbstractParameterEditor<?> editor = (AbstractParameterEditor<?>) param.getEditor();
				// update title font
				editor.setTitleFont(new Font(font.getFamily(), Font.ITALIC, font.getSize()));
				if (editor.getWidget() instanceof NumericTextField) {
					NumericTextField textField = (NumericTextField)editor.getWidget();
					textField.setFont(new Font(font.getFamily(), Font.ITALIC, font.getSize()));
					textField.setBackground(new Color(237,237,237));
					textField.setForeground(Color.DARK_GRAY);
				}
				editor.repaint();
			}	
		}
	}
	
	
	/* 
	 * Update the map parameter panel to reflect current map options
	 */
	private void updateMapPanel() {
		setChangeListenerEnabled(false);
		if (mapTypeParam.getValue() == MapType.PROBABILITIES) {
			setEnabledStyle(mapPOEParam, false);
			setEnabledStyle(mapLevelParam, true);
		
			if (intensityTypeParam.getValue() == IntensityType.PGA) {
				mapLevelParam.setConstraint(new DoubleConstraint(3,100));
				mapLevelParam.setValue(10);
				AbstractParameterEditor<?> editor = (AbstractParameterEditor<?>) mapLevelParam.getEditor();
				editor.setTitle("%g");
				editor.repaint();
				editor.refreshParamEditor();
			} else if (intensityTypeParam.getValue() == IntensityType.PGV) {
				mapLevelParam.setConstraint(new DoubleConstraint(2,100));
				mapLevelParam.setValue(10);
				AbstractParameterEditor<?> editor = (AbstractParameterEditor<?>) mapLevelParam.getEditor();
				editor.setTitle("cm/s");
				editor.repaint();
				editor.refreshParamEditor();
			} else if (intensityTypeParam.getValue() == IntensityType.PSA) {
				mapLevelParam.setConstraint(new DoubleConstraint(3,100));
				mapLevelParam.setValue(10);
				AbstractParameterEditor<?> editor = (AbstractParameterEditor<?>) mapLevelParam.getEditor();
				editor.setTitle("%g");
				editor.repaint();
				editor.refreshParamEditor();
			} else if (intensityTypeParam.getValue() == IntensityType.MMI) {
				mapLevelParam.setConstraint(new DoubleConstraint(4,10)); 
				AbstractParameterEditor<?> editor = (AbstractParameterEditor<?>) mapLevelParam.getEditor();
				mapLevelParam.setValue(6);
				editor.setTitle("MMI");
				editor.repaint();
				editor.refreshParamEditor();
			}
			
			mapLevelParam.getEditor().refreshParamEditor();
			mapPOEParam.getEditor().refreshParamEditor();
			
			
		} else if (mapTypeParam.getValue() == MapType.LEVEL) {
			setEnabledStyle(mapPOEParam, true);
			setEnabledStyle(mapLevelParam, false);
		}
		setChangeListenerEnabled(true);
	}
	
	
	/* 
	 * Utility functions for updating dynamic elements of the GUI, plot specifications, etc.
	 */
	private void updateRegionParamList(RegionType type, RegionCenterType centerType) {
		regionList.clear();
		regionList.addParameter(regionTypeParam);
		
		switch (type) {
		case CIRCULAR:
			regionList.addParameter(radiusParam);
			break;
		case CIRCULAR_WC94:
			// do nothing
			break;
		case RECTANGULAR:
			regionList.addParameter(minLatParam);
			regionList.addParameter(maxLatParam);
			regionList.addParameter(minLonParam);
			regionList.addParameter(maxLonParam);
			break;

		default:
			throw new IllegalStateException("Unknown region type: "+type);
		}
		
		if (type == RegionType.CIRCULAR) {
			regionList.addParameter(regionCenterTypeParam);
			
			if (centerType == RegionCenterType.SPECIFIED){
				regionList.addParameter(regionCenterLatParam);
				regionList.addParameter(regionCenterLonParam);
			}
		}
		
		regionList.addParameter(minDepthParam);
		regionList.addParameter(maxDepthParam);
		
		regionEditParam.getEditor().refreshParamEditor();
	} 
	
	// for plots: how to scale the symbols by size?
	private EvenlyDiscretizedFunc magSizeFunc;

	private EvenlyDiscretizedFunc getMagSizeFunc() {
			if (magSizeFunc != null)
				return magSizeFunc;
			
			// size function
			double minMag = 1.25;
			double magDelta = 0.5;
			int numMag = 2*8;
			magSizeFunc = new EvenlyDiscretizedFunc(minMag, numMag, magDelta);
	
			double dS = 3d;
			for (int i=0; i<magSizeFunc.size(); i++) {
				double mag = magSizeFunc.getX(i);
	
				// scale with stress drop, from Nicholas via e-mail 10/26/2015
	//			double radius = Math.pow((7d/16d)*Math.pow(10, 1.5*mag + 9)/(dS*1e6), 1d/3d) / 300d;
				// modified to have larger minimum size
				double radius = 2 + Math.pow((7d/16d)*Math.pow(10, 1.5*mag + 9)/(dS*1e6), 1d/3d) / 1000d;
				magSizeFunc.set(i, radius);
			}
			return magSizeFunc;
		}

	// for plotting: how to color the events by magnitude
	private CPT getMagCPT(){
		return getDistCPT();
	}
	
	// for plotting: how to color the events by distance
	private EvenlyDiscretizedFunc distFunc;

	private EvenlyDiscretizedFunc getDistFunc() {
		if (distFunc == null)
			distFunc = HistogramFunction.getEncompassingHistogram(0d, 366, 1d);
		return distFunc; 
	}

	private CPT distCPT;
	// for plotting: how to color events by distance
	private CPT getDistCPT() {
	
		if (distCPT != null){
			return distCPT;
		}
		
		EvenlyDiscretizedFunc distFunc = getDistFunc();
		double halfDelta = 0.5*distFunc.getDelta();
		
		int[][] colorMapValues = new int[][]{
			{68, 1, 84},
			{71, 13, 96},
			{72, 24, 106},
			{72, 35, 116},
			{71, 45, 123},
			{69, 55, 129},
			{66, 64, 134},
			{62, 73, 137},
			{59, 82, 139},
			{55, 91, 141},
			{51, 99, 141},
			{47, 107, 142},
			{44, 114, 142},
			{41, 122, 142},
			{38, 130, 142},
			{35, 137, 142},
			{33, 145, 140},
			{31, 152, 139},
			{31, 160, 136},
			{34, 167, 133},
			{40, 174, 128},
			{50, 182, 122},
			{63, 188, 115},
			{78, 195, 107},
			{94, 201, 98},
			{112, 207, 87},
			{132, 212, 75},
			{152, 216, 62},
			{173, 220, 48},
			{194, 223, 35},
			{216, 226, 25},
			{236, 229, 27},
			{251, 231, 35}
		};
		
		Color[] colorMap = new Color[colorMapValues.length];
		
		for (int i = 0; i < colorMap.length; i++){
			colorMap[i] = new Color(colorMapValues[i][0],colorMapValues[i][1],colorMapValues[i][2]);
		}
		
		distCPT = new CPT(distFunc.getMinX()-halfDelta, distFunc.getMaxX()+halfDelta, colorMap);
		return distCPT;
	}

//	private CPT probCPT;
	// for plotting: how to color events by probability
	private CPT getProbCPT() {
		CPT probCPT;
		
		// try these:
//		int[][] colorMapValues = new int[][]{
//			{115, 210, 230},
//			{109, 196, 227},
//			{104, 183, 224},
//			{98, 169, 222},
//			{92, 155, 219},
//			{87, 141, 217},
//			{82, 125, 214},
//			{77, 109, 212},
//			{71, 93, 209},
//			{66, 77, 207},
//			{61, 61, 204},
//			{75, 58, 207},
//			{89, 54, 209},
//			{103, 50, 212},
//			{116, 47, 214},
//			{130, 43, 217},
//			{150, 39, 212},
//			{170, 35, 208},
//			{190, 31, 204},
//			{210, 27, 199},
//			{230, 23, 195},
//		};
		

		int[][] colorMapValues = new int[][]{
		    {26,   212,   168},
		    {19,   203,   190},
		    {13,   194,   212},
		     {9,   188,   226},
		     {4,   182,   241},
		     {0,   176,   255},
		    {18,   165,   250},
		    {36,   154,   244},
		    {55,   143,   239},
		    {73,   132,   233},
		    {91,   121,   228},
		   {109,   110,   222},
		   {128,    99,   217},
		   {143,    90,   212},
		   {159,    80,   207},
		   {175,    71,   202},
		   {191,    61,   198},
		   {207,    52,   193},
		   {223,    42,   188},
		   {239,    33,   183},
		   {255,    23,   178},
		};
		
		Color[] colorMap = new Color[colorMapValues.length];
		
		for (int i = 0; i < colorMap.length; i++){
			colorMap[i] = new Color(colorMapValues[i][0],colorMapValues[i][1],colorMapValues[i][2]);
		}
		
		probCPT = new CPT(0, 1, colorMap);
			
		probCPT.setBlender(new Blender(){
			public Color blend(Color minColor, Color maxColor, float bias){
				return minColor;
			}
		});
		
		return probCPT;
	}
	
	private void buildFuncsCharsForBinned(XY_DataSet[] binnedFuncs,
			List<PlotElement> funcs, List<PlotCurveCharacterstics> chars, PlotSymbol sym) {
		EvenlyDiscretizedFunc magSizeFunc = getMagSizeFunc();
		double magDelta = magSizeFunc.getDelta();
		CPT magColorCPT = getMagCPT();
		
		for (int i=0; i<binnedFuncs.length; i++) {
			double mag = magSizeFunc.getX(i);
			XY_DataSet xy = binnedFuncs[i];
			if (xy.size() == 0)
				continue;
			xy.setName((float)(mag-0.5*magDelta)+" < M < "+(float)(mag+0.5*magDelta)
					+": "+xy.size()+" EQ");
			float size = (float)magSizeFunc.getY(i);
			Color c = magColorCPT.getColor((float)magSizeFunc.getX(i));
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(sym, size, c));
		}
	}

	/**
	 * Uses scalars from func2 to color funcs
	 * 
	 * @param binnedFuncs
	 * @param funcs
	 * @param chars
	 * @param cpt
	 * @param name2
	 * @param func2
	 */
	private void buildFuncsCharsForBinned2D(XY_DataSet[][] binnedFuncs, List<PlotElement> funcs,
			List<PlotCurveCharacterstics> chars, CPT cpt, String name2, EvenlyDiscretizedFunc func2, PlotSymbol sym) {
		EvenlyDiscretizedFunc magSizeFunc = getMagSizeFunc();
		double magDelta = magSizeFunc.getDelta();
		double func2Delta = func2.getDelta();
		for (int i=0; i<binnedFuncs.length; i++) {
			double mag = magSizeFunc.getX(i);
			for (int j=0; j<binnedFuncs[i].length; j++) {
				XY_DataSet xy = binnedFuncs[i][j];
				if (xy.size() == 0)
					continue;
				double scalar2 = func2.getX(j);
				String name = (float)(mag-0.5*magDelta)+" < M < "+(float)(mag+0.5*magDelta);
				name += ", "+(float)(scalar2-0.5*func2Delta)+" < "+name2+" < "+(float)(scalar2+0.5*func2Delta);
				name += ": "+xy.size()+" EQ";
				xy.setName(name);
				float size = (float)magSizeFunc.getY(i);
				Color c = cpt.getColor((float)func2.getX(j));
				funcs.add(xy);
				chars.add(new PlotCurveCharacterstics(sym, size, c));
			}
		}
	}

	private static final int tickLabelFontSize = 16;
	private static final int axisLabelFontSize = 16;
	private static final int plotLabelFontSize = 16;

	private static Color[] extra_colors = {Color.GRAY, Color.BLUE, Color.ORANGE, Color.GREEN};

	private static SimpleDateFormat catDateFormat = new SimpleDateFormat("yyyy\tMM\tdd\tHH\tmm\tss");
	private static final TimeZone utc = TimeZone.getTimeZone("UTC");
	static {
		catDateFormat.setTimeZone(utc);
	}

	private static void setupGP(GraphWidget widget) {
		widget.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, tickLabelFontSize));
		widget.setPlotLabelFontSize(plotLabelFontSize);
		widget.setAxisLabelFontSize(axisLabelFontSize);
		widget.setTickLabelFontSize(tickLabelFontSize);
		widget.setBackgroundColor(Color.WHITE);
	}
	// END Utility functions for updating dynamic elements of the GUI, plot specifications, etc.
	
	/*
	 * Begin functions for really doing stuff!
	 */
	private void fetchEvents() {
	
		
		String eventID = eventIDParam.getValue();

		if (eventID == null || eventID.isEmpty())
			System.err.println("Must supply event ID");
		else {
			mainshock = accessor.fetchEvent(eventID);
			if (mainshock == null){
				System.err.println("Event not found: " + eventID);
				setEnableParamsPostFetch(false);
			} else {
//				//fetch details
				ListIterator<?> iter = mainshock.getAddedParametersIterator();
				if(verbose){
					if (iter != null)
						while (iter.hasNext()){
							Parameter<?> param = (Parameter<?>) iter.next();
							System.out.println(param.getName() +":"+ param.getValue());
						}
				}
					
				System.out.println("Mainshock Mag/Lat/Lon/Depth: " + mainshock.getMag() + " " + mainshock.getHypocenterLocation());

				Double minDepth = minDepthParam.getValue();
				validateParameter(minDepth, "min depth");
				Double maxDepth = maxDepthParam.getValue();
				validateParameter(maxDepth, "max depth");

				// populate/validate data and forecast time windows
				updateForecastTimes();
				double dataMinDays = dataStartTimeParam.getValue();
				double dataMaxDays = dataEndTimeParam.getValue();
				
				
				if (regionTypeParam.getValue().isCircular()
						&& regionCenterTypeParam.getValue() == RegionCenterType.CENTROID) {

					// get initial region
					region = buildRegion(mainshock, mainshock.getHypocenterLocation());
					aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, region, region.getPlotWrap());

					ObsEqkRupList bigAftershocks = aftershocks.getRupsAboveMag(mainshock.getMag());
					largestShock = mainshock;
					while (!bigAftershocks.isEmpty()){
						System.out.println("Found an aftershock larger than the mainshock. Expanding aftershock zone...");
						// find the largest shock in the list
						double maxMag = largestShock.getMag();
						int maxMagIndex = 0;
						for(int i = 0 ; i < bigAftershocks.size(); i++){
							if(bigAftershocks.get(i).getMag() > maxMag){
								maxMag = bigAftershocks.get(i).getMag();
								maxMagIndex = i; 
							}
						}
						largestShock = bigAftershocks.get(maxMagIndex);

						// update region around largest shock
						region = buildRegion(largestShock, largestShock.getHypocenterLocation());
						aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, region, region.getPlotWrap());

						// look again for even larger shocks
						bigAftershocks = aftershocks.getRupsAboveMag(largestShock.getMag() + 0.1);
					}

					// now find centroid 
					if (aftershocks.isEmpty()) {
						if(verbose) System.out.println("No aftershocks found, skipping centroid...");
					} else {
						region = buildRegion(largestShock, getCentroid());
						aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, region, region.getPlotWrap());
					}

				} else {
					region = buildRegion(mainshock, null);
					aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, region, region.getPlotWrap());
				}
				
				// limit the catalog to MAX_EARTHQUAKE_NUMBER by changing the mc
				if (aftershocks.size() > MAX_EARTHQUAKE_NUMBER){
					System.out.println("Found " + aftershocks.size() + " aftershocks. " + 
							" Keeping only the " + MAX_EARTHQUAKE_NUMBER + " largest.");
					double[] magnitudes = ETAS_StatsCalc.getAftershockMags(aftershocks);
					Arrays.sort(magnitudes);
					double newMc = magnitudes[aftershocks.size() - MAX_EARTHQUAKE_NUMBER];
					aftershocks = aftershocks.getRupsAboveMag(newMc);
				}
				
				System.out.println("Found " + aftershocks.size() + " aftershocks after filtering.");

				//update the Aftershock Zone editor with the collection radius and other info
				double newRadius = (region.getMaxLat() - region.getMinLat())*111.1111/2;
				if(D) System.out.println(region + "\n\tRadius: " + newRadius);
				
				radiusParam.setValue(newRadius);
				radiusParam.getEditor().refreshParamEditor();
				minLatParam.setValue(region.getMinLat());
				minLatParam.getEditor().refreshParamEditor();
				maxLatParam.setValue(region.getMaxLat());
				maxLatParam.getEditor().refreshParamEditor();
				minLonParam.setValue(region.getMinLon());
				minLonParam.getEditor().refreshParamEditor();
				maxLonParam.setValue(region.getMaxLon());
				maxLonParam.getEditor().refreshParamEditor();
				
				regionCenterLatParam.setValue((region.getMaxLat()+region.getMinLat())/2d);
				regionCenterLatParam.getEditor().refreshParamEditor();
				regionCenterLonParam.setValue((region.getMaxLon()+region.getMinLon())/2d);
				regionCenterLonParam.getEditor().refreshParamEditor();
			}
		}
	}
	
	private SphRegion buildRegion(ObsEqkRupture event, Location centroid) {
		SphRegion result;
		RegionType type = regionTypeParam.getValue();
		
		if (type == RegionType.CIRCULAR || type == RegionType.CIRCULAR_WC94) {
			double radius;
			if (type == RegionType.CIRCULAR) {
				radius = radiusParam.getValue();
			} else {
				if (wcMagLen == null)
					wcMagLen = new WC1994_MagLengthRelationship();
				radius = 10d + 1.5*wcMagLen.getMedianLength(event.getMag()); //multiplied by 1.5 because length is already 2x radius
				if(verbose) System.out.println("Collecting Aftershocks within 10 km + 3 WC94 Radii: "+ (float)radius + " km");
			}
			
			RegionCenterType centerType = regionCenterTypeParam.getValue();
			Location loc;
			if (centerType == RegionCenterType.EPICENTER)
				loc = event.getHypocenterLocation();
			else if (centerType == RegionCenterType.CENTROID)
				loc = centroid;
			else if (centerType == RegionCenterType.SPECIFIED)
				loc = new Location(regionCenterLatParam.getValue(), regionCenterLonParam.getValue());
			else
				throw new IllegalStateException("Unknown Region Center Type: "+centerType);
			
			//return new Region(loc, radius);
			result = SphRegion.makeCircle (new SphLatLon(loc), radius);
		} else  if (type == RegionType.RECTANGULAR) {
			Location lower = new Location(minLatParam.getValue(), minLonParam.getValue());
			Location upper = new Location(maxLatParam.getValue(), maxLonParam.getValue());
			//return new Region(lower, upper);
			result = SphRegion.makeMercRectangle (new SphLatLon(lower), new SphLatLon(upper));
		} else {
			throw new IllegalStateException("Unknown region type: "+type);
		}

		// If the event (i.e. mainshock) is outside the plotting domain, change its
		// hypocenter so it is inside the plotting domain

		Location hypo = event.getHypocenterLocation();
		if (result.getPlotWrap()) {
			if (hypo.getLongitude() < 0.0) {
				event.setHypocenterLocation (new Location (
					hypo.getLatitude(), hypo.getLongitude() + 360.0, hypo.getDepth() ));
			}
		} else {
			if (hypo.getLongitude() > 180.0) {
				event.setHypocenterLocation (new Location (
					hypo.getLatitude(), hypo.getLongitude() - 360.0, hypo.getDepth() ));
			}
		}

		return result;
	}
	
	private void updateForecastTimes(){
		double elapsedDays = (double) (System.currentTimeMillis() - mainshock.getOriginTime())/ETAS_StatsCalc.MILLISEC_PER_DAY;

		double forecastStartDays;
		Double dataMaxDays;
		if (nowBoolean.getValue() || (forecastStartTimeParam.getValue() == null)){
			if (elapsedDays > dataEndTimeParam.getMax()){
				dataMaxDays = dataEndTimeParam.getMax();
				forecastStartDays = dataEndTimeParam.getMax();
			} else {
				dataMaxDays = elapsedDays;
				forecastStartDays = elapsedDays;
			}

		} else if (dataEndTimeParam.getValue() == null){
			forecastStartDays = forecastStartTimeParam.getValue();
			dataMaxDays = forecastStartDays;

		} else {
			// set up to link endtime and start time.
			forecastStartDays = forecastStartTimeParam.getValue();
			dataMaxDays = forecastStartDays;
			
		}

		forecastStartTimeParam.setValue(round(forecastStartDays));
		forecastStartTimeParam.getEditor().refreshParamEditor();

		forecastEndTimeParam.setValue(round(forecastStartDays + forecastDurationParam.getValue()));
		forecastEndTimeParam.getEditor().refreshParamEditor();

		//check for data end time greater than current time and set to current time if smaller
		if (elapsedDays < dataMaxDays){
			if(verbose) System.out.println("Setting data end time and forecast start time to current time.");
			dataMaxDays = elapsedDays;
		}
		validateParameter(dataMaxDays, "end time");
		validateParameter(forecastStartDays, "forecast start time");

		dataEndTimeParam.setValue(round(dataMaxDays));
		dataEndTimeParam.getEditor().refreshParamEditor();				

		if(verbose) System.out.println("dataStartTime = " + dataStartTimeParam.getValue() + 
				" dataEndTime = " + dataEndTimeParam.getValue() +
				" elapsedTime = " + elapsedDays);
		if(verbose) System.out.println("forecastStart = " + forecastStartTimeParam.getValue() +
				" forecastEnd = " + forecastEndTimeParam.getValue() +
				" forecastDuration = " + forecastDurationParam.getValue());
	}
	
	private void wrapLongitudes(boolean wrap){
		double lon;
		
		if (wrap){
			try{
				lon = regionCenterLonParam.getValue();
				if(lon < 0){
					lon+=360;
					regionCenterLonParam.setValue(lon);
					regionCenterLonParam.getEditor().refreshParamEditor();
				}
			} catch (Exception e){
				//doesn't exist
			}
			
			try{
				lon = minLonParam.getValue();
				if(lon < 0){
					lon+=360;
					minLonParam.setValue(lon);
					minLonParam.getEditor().refreshParamEditor();
				}
			} catch (Exception e) {
				//doesn't exist
			}

			try{
				lon = maxLonParam.getValue();
				if(lon < 0){
					lon+=360;
					maxLonParam.setValue(lon);
					maxLonParam.getEditor().refreshParamEditor();
				}
			} catch (Exception e) {
				//doesn't exist
			}
			
		}
	}
	
	private void setMainshock(ObsEqkRupture mainshock) {
		if(verbose) System.out.println("Setting up Mainshock...");
		if (mainshock == null){
			return;
		}
		this.mainshock = mainshock;
		
		// if the mainshock has a longitude > 180, change the longitude parameters to match the wrapping
		wrapLongitudes(mainshock.getHypocenterLocation().getLongitude() > 180);
		
		genericParams = null;
		tectonicRegimeParam.setValueAsDefault();
		
		TectonicRegime regime;
		
		try {
			if (genericFetch == null)
				genericFetch = new GenericETAS_ParametersFetch();
			
			if(verbose) System.out.println("Determining tectonic regime...");
			consoleScroll.repaint();

			TectonicRegimeTable regime_table = new TectonicRegimeTable();
			String regimeName = regime_table.get_strec_name(mainshock.getHypocenterLocation().getLatitude(), mainshock.getHypocenterLocation().getLongitude());
			regime = TectonicRegime.forName(regimeName);
					
		} catch (RuntimeException e) {
			System.out.println("Error fetching generic params. Assigning global average values");
			regime = TectonicRegime.GLOBAL_AVERAGE;	
		}
		
		if (regime == null){
			System.out.println("Error fetching generic params. Assigning global average values");
			regime = TectonicRegime.GLOBAL_AVERAGE;
		}			
				
		genericParams = genericFetch.get(regime);
		try{
			Preconditions.checkNotNull(genericParams, "Generic params not found or server error");
		}catch (Exception e){
			System.out.println("Error retrieving generic parameters for tectonic regime: " + regime);
		}
		if(verbose) System.out.println("Generic params for "+regime+": "+genericParams);
		
		aValRangeParam.setValue(new Range(round(genericParams.get_a()-3*genericParams.get_aSigma(),2), round(genericParams.get_a()+3*genericParams.get_aSigma(),2)));
		pValRangeParam.setValue(new Range(round(genericParams.get_p()-3*genericParams.get_pSigma(),2), round(genericParams.get_p()+3*genericParams.get_pSigma(),2)));
		cValRangeParam.setValue(new Range(round(genericParams.get_c()/Math.pow(10, 3*genericParams.get_logcSigma()),sigDigits), Math.min(1, round(genericParams.get_c()*Math.pow(10, 3*genericParams.get_logcSigma()),sigDigits))));
		
		bParam.setValue(round(genericParams.get_b(),2));
		bParam.getEditor().refreshParamEditor();
		
		magRefParam.setValue(round(genericParams.get_refMag(),2));
		magRefParam.getEditor().refreshParamEditor();
		
		tectonicRegimeParam.setValue(regime);
		tectonicRegimeParam.getEditor().refreshParamEditor();
		
		setMagComplete();
		resetFitConstraints(genericParams);
		updateForecastTimes();
		
		// as a courtesy, spit out the decimal days remaining in the origin day
		System.out.println("The mainshock occurred " + String.format("%.4f", getTimeRemainingInUTCDay()) + " days before midnight (UTC)\n");
		
		if(verbose) System.out.println("Building generic model...");
		// initialize the generic model using an assumed global Mmax of 9.5, max generation depth 100, number of simulations 0
		genericModel = new ETAS_AftershockModel_Generic(mainshock, aftershocks, genericParams, 
				dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
				forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), mcParam.getValue(),
				9.5, 100, 0, fitMSProductivityParam.getValue(), timeDepMcParam.getValue());

	} //end setupGUI
	
	private void setMagComplete(double mc){
		double b_value = 1; 
		if (ogataMND != null){
			 b_value = ogataMND.calculateBWithFixedMc(mc);
			 if(D) System.out.println("mc: " + mc + " b: " + b_value);
		}

		mcParam.setValue(round(mc,2));
		mcParam.getEditor().refreshParamEditor();
		bParam.setValue(round(b_value,2));
		bParam.getEditor().refreshParamEditor();

		link_alpha(); //reset alpha based on b
		
//		printTip(2);

	}
	
	private void setMagComplete(){
		// this routine fits by b-value stability. Ogata&Katsura simultaneous Mc and b is set with a different routine.
		
		aftershockMND = ObsEqkRupListCalc.getMagNumDist(aftershocks, -3.05, 131, 0.1);
		double mc = 0;
		double b_value = 1.0;
		double b_sigma = 0.1;
		
		// do it by the Ogata and Katsura method
		if (aftershocks.size() > 0){
			if (genericParams != null){
				b_value = genericParams.get_b();
				b_sigma = genericParams.get_bSigma();
			}
			ogataMND = new OgataMagFreqDist(aftershocks, b_value, b_sigma); //compute entire-magnitude-range Mc with b fixed 
			
			// iterate to a stable estimate solution
			mc = ogataMND.calculateMcWithFixedB(b_value);
			if(D) System.out.println("mc: " + mc + " b: " + b_value);

		} else {
		 	mc = genericParams.get_refMag() + 0.5;
		}
		
		mcParam.setValue(round(mc,2));
		mcParam.getEditor().refreshParamEditor();
		bParam.setValue(round(b_value,2));
		bParam.getEditor().refreshParamEditor();
		
		link_alpha(); //reset alpha based on b
		
//		printTip(2);
	}
	
	private void plotAftershockHypocs() {
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		EvenlyDiscretizedFunc magSizeFunc = getMagSizeFunc();
		
		boolean colorByTime = true;
		CPT timeCPT = null;
		
		if (colorByTime) {
			timeCPT = getDistCPT().rescale(0d, dataEndTimeParam.getValue() + 1e-6);
		}
		PaintScaleLegend subtitle = null;
		
		RuptureSurface mainSurf = mainshock.getRuptureSurface();
		if (mainSurf != null && !mainSurf.isPointSurface()) {
			FaultTrace trace = mainshock.getRuptureSurface().getEvenlyDiscritizedUpperEdge();
			DefaultXY_DataSet traceFunc = new DefaultXY_DataSet();
			traceFunc.setName("Main Shock Trace");
			for(Location loc:trace){
				traceFunc.set(loc.getLongitude(), loc.getLatitude());
			}
			funcs.add(traceFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		} else {
			Location hypo = mainshock.getHypocenterLocation();
			DefaultXY_DataSet xy = new DefaultXY_DataSet(new double[] {hypo.getLongitude()},
					new double[] {hypo.getLatitude()});
			xy.setName("Main Shock Location");
			funcs.add(xy);
			float size = (float)magSizeFunc.getY(magSizeFunc.getClosestXIndex(mainshock.getMag()));
			Color c;
			if (colorByTime)
				c = timeCPT.getMinColor();
			else
				c = Color.BLACK;
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, size, c));
		}
		
		// now aftershocks
		if(aftershocks.size() > 0){
			List<Point2D> points = Lists.newArrayList();
			List<Double> mags = Lists.newArrayList();
			List<Double> timeDeltas = Lists.newArrayList();
			for (ObsEqkRupture rup : aftershocks) {
				Location loc = rup.getHypocenterLocation();
				points.add(new Point2D.Double(loc.getLongitude(), loc.getLatitude()));
				mags.add(rup.getMag());
				timeDeltas.add(getTimeSinceMainshock(rup));
			}

			if (colorByTime) {
				EvenlyDiscretizedFunc timeFunc = HistogramFunction.getEncompassingHistogram(0d, timeCPT.getMaxValue()*0.99, 0.1d);
				XY_DataSet[][] aftershockDatasets = XY_DatasetBinner.bin2D(points, mags, timeDeltas, magSizeFunc, timeFunc);

				buildFuncsCharsForBinned2D(aftershockDatasets, funcs, chars, timeCPT, "time", timeFunc, PlotSymbol.CIRCLE);

				double cptInc = 0d;
				if ((timeCPT.getMaxValue() - timeCPT.getMinValue()) < 10)
					cptInc = 1d;
				subtitle = XYZGraphPanel.getLegendForCPT(timeCPT, "Time (days)", axisLabelFontSize, tickLabelFontSize,
						cptInc, RectangleEdge.RIGHT);
			} else {
				XY_DataSet[] aftershockDatasets = XY_DatasetBinner.bin(points, mags, magSizeFunc);

				buildFuncsCharsForBinned(aftershockDatasets, funcs, chars, PlotSymbol.CIRCLE);
			}
		}
		
		// now add outline
		if (region != null) {
			DefaultXY_DataSet outline = new DefaultXY_DataSet();
			for (Location loc : region.getBorder())
				outline.set(loc.getLongitude(), loc.getLatitude());
			outline.set(outline.get(0));
			outline.setName("Region Outline");
			
			funcs.add(outline);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		}
		
		Collections.reverse(funcs);
		Collections.reverse(chars);
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Aftershock Epicenters", "Longitude", "Latitude");
				
 		if (subtitle != null)
  			spec.addSubtitle(subtitle);
 		
		if (epicenterGraph == null){
			epicenterGraph = new GraphWidget(spec);
			((NumberAxis) epicenterGraph.getGraphPanel().getXAxis()).setAutoRangeIncludesZero(false);;
			((NumberAxis) epicenterGraph.getGraphPanel().getYAxis()).setAutoRangeIncludesZero(false);;
			
			Component buttonPanel = epicenterGraph.getButtonControlPanel();
			buttonPanel.setVisible(false);
			GraphPanel graphPanel = epicenterGraph.getGraphPanel();
			graphPanel.getComponent(2).setVisible(false);
		} else
			epicenterGraph.setPlotSpec(spec);	
		
		double regBuff = 0.05;
		if (region != null) {
			epicenterGraph.setAxisRange(region.getMinLon()-regBuff, region.getMaxLon()+regBuff,
					region.getMinLat()-regBuff, region.getMaxLat()+regBuff);
		} else {
			 MinMaxAveTracker latTrack = new MinMaxAveTracker();
			 MinMaxAveTracker lonTrack = new MinMaxAveTracker();
			 latTrack.addValue(mainshock.getHypocenterLocation().getLatitude());
			 lonTrack.addValue(mainshock.getHypocenterLocation().getLongitude());
			 for (ObsEqkRupture rup : aftershocks) {
				 Location loc = rup.getHypocenterLocation();
				 latTrack.addValue(loc.getLatitude());
				 lonTrack.addValue(loc.getLongitude());
			 }
			 epicenterGraph.setAxisRange(lonTrack.getMin()-regBuff, lonTrack.getMax()+regBuff,
						latTrack.getMin()-regBuff, latTrack.getMax()+regBuff);
		}
		
		setupGP(epicenterGraph);
 

      			
			
		
	}
	
	private void plotMFDs(IncrementalMagFreqDist mfd, double magComplete) {
		if(D) System.out.println("start plotMFDs...");
		List<XY_DataSet> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		 
		double plotMinMag = Double.POSITIVE_INFINITY;
		double plotMaxMag = mainshock.getMag();
		for (int i=0; i<mfd.size(); i++) {
			if (mfd.getY(i) > 0) {
				double mag = mfd.getX(i);
				plotMinMag = Math.min(plotMinMag, mag);
				plotMaxMag = Math.max(plotMaxMag, mag);
			}
		}
		if (Double.isInfinite(plotMinMag))
			plotMinMag = 0d;
		plotMinMag = Math.floor(plotMinMag);
		plotMaxMag = Math.ceil(plotMaxMag);
		
		EvenlyDiscretizedFunc cmlMFD =  mfd.getCumRateDistWithOffset();
		cmlMFD.setName("Observed");
		if (aftershocks.size() > 0){
			funcs.add(cmlMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		}
		
		ArbitrarilyDiscretizedFunc cmlPoints = new ArbitrarilyDiscretizedFunc();
		cmlPoints.setName("");

		double prevVal = cmlMFD.getY(0);

		for (int i=1; i<cmlMFD.size(); i++) {
			double val = cmlMFD.getY(i);
			if (val != prevVal) {
				cmlPoints.set(cmlMFD.getX(i-1), prevVal);

				prevVal = val;
			}
		}
		if (aftershocks.size() > 0){
			funcs.add(cmlPoints);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 6f, Color.BLACK));
		}
		
		// add Ogata entire-magnitude-range fit
			if (ogataMND != null){
				EvenlyDiscretizedFunc fitEMR = ogataMND.getCumulativeMND(plotMinMag, genericParams.maxMag, magPrecisionParam.getValue());
				fitEMR.setName("Model"); 
				funcs.add(fitEMR);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, bayesian_color));
			}
		//
			
		double plotMinY = 0.9d;
		double plotMaxY = Math.max(10d,(cmlMFD.getMaxY() + 1d)*10d);
		
		List<Double> yValsForVerticalLines = Lists.newArrayList(0d, 1e-16, plotMinY, 1d, plotMaxY, 1e3, 2e3, 3e3 ,4e3, 5e3);
		Collections.sort(yValsForVerticalLines);
		
		// add mainshock mag
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		for (double y : yValsForVerticalLines)
			xy.set(mainshock.getMag(), y);
		xy.setName(String.format("Mainshock (M%2.1f)", mainshock.getMag()));
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, generic_color));
		
		if(verbose){ System.out.println("Calculating MFD with b: "+bParam.getValue()); }
		
		if (bParam.getValue() != null) {
			// add Mc used for b-value calculation
			double mc = mcParam.getValue();
			xy = new DefaultXY_DataSet();
			for (double y : yValsForVerticalLines)
				xy.set(mc, y);
			xy.setName(String.format("Mag complete (M%2.1f)", mc));
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, sequence_specific_color));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Magnitude Distribution", "Magnitude", "Count");
		spec.setLegendVisible(true);
		
		if(D) System.out.println("plot MFD: " + spec.getTitle());
		
		
		if (magNumGraph == null){
			magNumGraph = new GraphWidget(spec);
			Component buttonPanel = magNumGraph.getButtonControlPanel();
			buttonPanel.setVisible(true);
			GraphPanel graphPanel = magNumGraph.getGraphPanel();
			graphPanel.getComponent(2).setVisible(false);
			
			
		} else
			magNumGraph.setPlotSpec(spec);

		magNumGraph.setY_Log(true);

		magNumGraph.setY_AxisRange(plotMinY, plotMaxY);
		magNumGraph.setX_AxisRange(plotMinMag, plotMaxMag);
		setupGP(magNumGraph);

		if(D) System.out.println("...end plotMFDs");
	}
	
	private void plotMagVsTime() {
		List<Point2D> points = Lists.newArrayList();
		List<Double> mags = Lists.newArrayList();
		
		points.add(new Point2D.Double(getTimeSinceMainshock(mainshock), mainshock.getMag()));
		mags.add(mainshock.getMag());
		
		MinMaxAveTracker magTrack = new MinMaxAveTracker();
		magTrack.addValue(mainshock.getMag());
		for (int i=0; i<aftershocks.size(); i++) {
			ObsEqkRupture aftershock = aftershocks.get(i);
			points.add(new Point2D.Double(getTimeSinceMainshock(aftershock), aftershock.getMag()));
			mags.add(aftershock.getMag());
			magTrack.addValue(aftershock.getMag());
		}
		
		boolean colorByDist = true;
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		EvenlyDiscretizedFunc magSizeFunc = getMagSizeFunc();
		PaintScaleLegend subtitle = null;
		
		if (colorByDist) {
			List<Double> dists = Lists.newArrayList();
			// TODO horizontal, correct?
			dists.add(0d); // mainshock distance from itself, thus zero
			
			double maxDist = 1e-6;
			double thisDist;
			for (int i=0; i<aftershocks.size(); i++) {
				ObsEqkRupture aftershock = aftershocks.get(i);
				thisDist = LocationUtils.horzDistanceFast(mainshock.getHypocenterLocation(), aftershock.getHypocenterLocation());
				if (thisDist > maxDist){
					maxDist = thisDist;
				}
				dists.add(thisDist);
			}
			
			EvenlyDiscretizedFunc distFunc = getDistFunc();
			
			XY_DataSet[][] binnedFuncs = XY_DatasetBinner.bin2D(points, mags, dists, magSizeFunc, distFunc);
			
			CPT distCPT = getDistCPT().rescale(0, maxDist);
			
			buildFuncsCharsForBinned2D(binnedFuncs, funcs, chars, distCPT, "dist", distFunc, PlotSymbol.FILLED_CIRCLE);
			
			subtitle = XYZGraphPanel.getLegendForCPT(distCPT, "Distance (km)", axisLabelFontSize, tickLabelFontSize,
					0d, RectangleEdge.RIGHT);
		} else {
			XY_DataSet[] magBinnedFuncs = XY_DatasetBinner.bin(points, mags, magSizeFunc);
			
			buildFuncsCharsForBinned(magBinnedFuncs, funcs, chars, PlotSymbol.CIRCLE);
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Magnitude Vs Time", "Days Since Mainshock", "Magnitude");
		if (subtitle != null)
			spec.addSubtitle(subtitle);
		
		if (magTimeGraph == null){
			magTimeGraph = new GraphWidget(spec);
			GraphPanel graphPanel = magTimeGraph.getGraphPanel();
			graphPanel.getComponent(2).setVisible(false);
		}
		
		else
			magTimeGraph.setPlotSpec(spec);
		
		magTimeGraph.setX_AxisRange(0.001, dataEndTimeParam.getValue()+0.75);
		magTimeGraph.setY_AxisRange(Math.max(0, magTrack.getMin()-1d), magTrack.getMax()+1d);
		setupGP(magTimeGraph);
		
	
	}
	
	private void plotCumulativeNum() {
		if(verbose) System.out.println("Plotting cumulative number with time...");
		double magMin;
		
		magMin = mcParam.getValue();
		
		ArbitrarilyDiscretizedFunc countFunc = new ArbitrarilyDiscretizedFunc();
		double count = 0;
		
		aftershocks.sortByOriginTime();
		for (int i=0; i<aftershocks.size(); i++) {
			ObsEqkRupture aftershock = aftershocks.get(i);
			if (aftershock.getMag() < magMin)
				continue;
			double time = getTimeSinceMainshock(aftershock);
			count++;
			countFunc.set(time, count);
		}
		countFunc.set(dataEndTimeParam.getValue(), count);
		countFunc.setName("Data");
		
		double maxY = count + 1;
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(countFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
//		if (seqSpecModel != null) {
//			if (seqSpecModel.get_nSims() != 0 && devMode) {
//				ArbitrarilyDiscretizedFunc expected = getModelCumNumWithLogTimePlot(seqSpecModel, magMin);
//
//				maxY = Math.max(count, expected.getMaxY());
//
//				funcs.add(expected);
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, sequence_specific_color));
//				expected.setName("Seq. specific");
//				
//				ArbitrarilyDiscretizedFunc lower = getFractileCumNumWithLogTimePlot(seqSpecModel, magMin, 0.025);
//
//				funcs.add(lower);
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, sequence_specific_color));
//				lower.setName("");
//
//				ArbitrarilyDiscretizedFunc upper = getFractileCumNumWithLogTimePlot(seqSpecModel, magMin, 0.975);
//				maxY = Math.max(maxY, upper.getMaxY());
//
//				funcs.add(upper);
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, sequence_specific_color));
//				upper.setName("");
//
//				UncertainArbDiscDataset uncertainFunc = new UncertainArbDiscDataset(expected, lower, upper);
//				funcs.add(uncertainFunc);
//				uncertainFunc.setName("Seq. specific range");
//				PlotLineType plt = PlotLineType.SHADED_UNCERTAIN;
//				Color seq_spec_trans_color = new Color(sequence_specific_color.getRed(), sequence_specific_color.getGreen(),sequence_specific_color.getBlue(), (int) (0.3*255) );
//				chars.add(new PlotCurveCharacterstics(plt, 1f, seq_spec_trans_color));
//			}
//
//			ArbitrarilyDiscretizedFunc fit = getModelFit(seqSpecModel, magMin);
//			funcs.add(fit);
//			if(devMode){
//				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, sequence_specific_color));
//				fit.setName("Seq. specific Fit");
//			} else {
//				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, bayesian_color));
//				fit.setName("");
//			}
//		} 
		
		if (genericModel != null && !plotSpecificOnlyParam.getValue()) {
			if (genericModel.nSims != 0) {
				ArbitrarilyDiscretizedFunc expected = getFractileCumNumWithLogTimePlot(genericModel, magMin, 0.5);
				maxY = Math.max(maxY, expected.getMaxY());

				if (bayesianModel == null || devMode){
					funcs.add(expected);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, generic_color));
					expected.setName("Typical sequence");
				}

				ArbitrarilyDiscretizedFunc lower = getFractileCumNumWithLogTimePlot(genericModel, magMin, 0.025);

				funcs.add(lower);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, generic_color));
				lower.setName("");

				ArbitrarilyDiscretizedFunc upper = getFractileCumNumWithLogTimePlot(genericModel, magMin, 0.975);
				maxY = Math.max(maxY, upper.getMaxY());

				funcs.add(upper);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, generic_color));
				upper.setName("");

				UncertainArbDiscDataset uncertainFunc = new UncertainArbDiscDataset(expected, lower, upper);
				funcs.add(uncertainFunc);
				uncertainFunc.setName("Typical range");
				PlotLineType plt = PlotLineType.SHADED_UNCERTAIN;
				Color generic_color_trans = new Color(generic_color.getRed(), generic_color.getGreen(),generic_color.getBlue(), (int) (0.3*255) );
				chars.add(new PlotCurveCharacterstics(plt, 1f, generic_color_trans));
			}
//			if(devMode){
//				ArbitrarilyDiscretizedFunc fit = getModelFit(genericModel, magMin);
//				funcs.add(fit);
//				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, generic_color));
//				fit.setName("Typical sequence");
//			}
		}
		
		if (bayesianModel != null) {
			if (bayesianModel.nSims != 0) {
				ArbitrarilyDiscretizedFunc expected = getModelCumNumWithLogTimePlot(bayesianModel, magMin);
				maxY = Math.max(maxY, expected.getMaxY());

				funcs.add(expected);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, bayesian_color));
				expected.setName("This sequence");
				
				ArbitrarilyDiscretizedFunc lower = getFractileCumNumWithLogTimePlot(bayesianModel, magMin, 0.025);

				funcs.add(lower);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, bayesian_color));
				lower.setName("");

				ArbitrarilyDiscretizedFunc upper = getFractileCumNumWithLogTimePlot(bayesianModel, magMin, 0.975);
				maxY = Math.max(maxY, upper.getMaxY());

				funcs.add(upper);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, bayesian_color));
				upper.setName("");

				UncertainArbDiscDataset uncertainFunc = new UncertainArbDiscDataset(expected, lower, upper);
				funcs.add(uncertainFunc);
				uncertainFunc.setName("Range for this sequence");
				PlotLineType plt = PlotLineType.SHADED_UNCERTAIN;
				Color bayesian_color_trans = new Color(bayesian_color.getRed(), bayesian_color.getGreen(),bayesian_color.getBlue(), (int) (0.3*255) );
				chars.add(new PlotCurveCharacterstics(plt, 1f, bayesian_color_trans));
			}
			ArbitrarilyDiscretizedFunc fit = getModelFit(bayesianModel, magMin);
			funcs.add(fit);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, bayesian_color));
			fit.setName("");
		}

		PlotSpec spec = new PlotSpec(funcs, chars, "Cumulative number of aftershocks above M"+String.format("%2.1f", magMin), "Days Since Mainshock",
				"Cumulative Number of Aftershocks");
		spec.setLegendVisible(true);
		
		if (cmlNumGraph == null){
			cmlNumGraph = new GraphWidget(spec);
			GraphPanel graphPanel = cmlNumGraph.getGraphPanel();
			graphPanel.getComponent(2).setVisible(false);
		}
		else
			cmlNumGraph.setPlotSpec(spec);
		
		setupGP(cmlNumGraph);
		cmlNumGraph.setY_AxisRange(1d, Math.max(maxY,1)*1.1);
		
		if(D) System.out.println("... end plotCumulativeNum");
	}
	
	private ArbitrarilyDiscretizedFunc getModelFit(ETAS_AftershockModel model, double magMin){
		double tMin = dataStartTimeParam.getValue();
		double tMax = dataEndTimeParam.getValue();
		if (tMax > tMin+1e-9) tMax += -1e-9;	//this is related to preventing log of zero
		
		if(D) System.out.println("Plotting fit for "+ tMin +" - "+ tMax + " days");
		int numPts = 100;
		return model.getExpectedNumEventsWithLogTime(magMin, tMin, tMax, numPts);
	}

	private ArbitrarilyDiscretizedFunc getModelCumNumWithLogTimePlot(ETAS_AftershockModel model, double magMin) {
		return getFractileCumNumWithLogTimePlot(model, magMin, 0.5);
	}

	private ArbitrarilyDiscretizedFunc getFractileCumNumWithLogTimePlot(ETAS_AftershockModel model, double magMin, double fractile) {
		double tMin = model.getForecastMinDays();
		double tMax = model.getForecastMaxDays();
		
		Preconditions.checkState(tMax >= tMin);
		
		int numPts = 100;
		return model.getFractileCumNumEventsWithLogTime(magMin, tMin, tMax, numPts, fractile);
	}
	
	private void loadCatalog(File catalogFile) throws IOException {
			List<String> lines = Files.readLines(catalogFile, Charset.defaultCharset());
			ObsEqkRupList myAftershocks = new ObsEqkRupList();
			ObsEqkRupture myMainshock = null;
			Double catalogMaxDays = null;
			String catalogEventID = "custom";
			
			for (int i=0; i<lines.size(); i++) {
				String line = lines.get(i).trim();
				if (line.startsWith("#")) {
					if (line.toLowerCase().startsWith("# main")
							&& i < lines.size()-1 && lines.get(i+1).startsWith("#")) {
					
						// main shock on next line, starting with a #
						String mainshockLine = lines.get(i+1).substring(1).trim();
						if(verbose) System.out.println("Detected mainshock in file: "+mainshockLine);
						try {
							myMainshock = fromCatalogLine(mainshockLine);
						} catch (Exception e) {
							System.err.println("Error loading mainshock");
						}
					} else if (line.toLowerCase().startsWith("# header:")){
						String headerLine = line.substring(1).trim();
						if(verbose) System.out.println(headerLine);
						try {
							catalogEventID = fromHeaderLine(headerLine, "eventID");
							catalogMaxDays = Double.parseDouble(fromHeaderLine(headerLine, "dataEndTime"));
						} catch (Exception e) {
							if(D) e.printStackTrace();
							System.err.println("Error loading header information");
						}
						
					}
						
					continue;
				}
				if (!line.isEmpty()) {
					try {
						myAftershocks.add(fromCatalogLine(line));
					} catch (ParseException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
				}
			}
			
			System.out.println("Loaded "+myAftershocks.size()+" aftershocks from file.");
			
			if (myMainshock == null) {
				// no mainshock detected, must load or use existing
				boolean prompt = true;     
				if (mainshock != null) {
					// ask the user if they want to overwrite the existing one in the app
					int ret = JOptionPane.showConfirmDialog(this, "A main shock has already been loaded."
							+ "\nDo you wish to specify your own custom main shock instead?",
							"Specify Custom Mainshock?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
					prompt = ret == JOptionPane.YES_OPTION;
				}
				StringParameter lineParam = new StringParameter("Mainshock Line");
				lineParam.setValue("Year Month Day Hour Minute Sec Lat Lon Depth Magnitude");
				while (prompt) {
					try {
						myMainshock = fromCatalogLine(lineParam.getValue());
						break;
					} catch (Exception e) {
						int ret = JOptionPane.showConfirmDialog(this, "Error: "+e.getMessage()+"\nTry again?",
								"Error Parsing Mainshock", JOptionPane.OK_CANCEL_OPTION);
						if (ret == JOptionPane.CANCEL_OPTION)
							break;
					}
				}
			}
			
			Preconditions.checkState(myMainshock != null, "Could not laod mainshock");
			
			System.out.println("Mainshock Mag/Lat/Lon/Depth: " + myMainshock.getMag() + " " + myMainshock.getHypocenterLocation());
			
			if (catalogMaxDays == null){
				if (verbose) System.out.println("No catalog time found in header. Using time of last aftershock.");
				if (myAftershocks.size() > 0)
					catalogMaxDays = (double) (myAftershocks.get(myAftershocks.size()-1).getOriginTime() - myMainshock.getOriginTime() )/ETAS_StatsCalc.MILLISEC_PER_DAY;
				else
					catalogMaxDays = 0d;
			}
			
			if ( dataEndTimeParam.getValue() == null || catalogMaxDays < dataEndTimeParam.getValue() ) {
				dataEndTimeParam.setValue(catalogMaxDays);
				dataEndTimeParam.getEditor().refreshParamEditor();
			}
			
			mainshock = myMainshock;
			aftershocks = myAftershocks.getRupsBefore((long) (mainshock.getOriginTime() + dataEndTimeParam.getValue()*ETAS_StatsCalc.MILLISEC_PER_DAY));
			
			//rebuild a region based on the largest event
			largestShock = mainshock;
			double maxMag = mainshock.getMag();
			for(int i = 0 ; i < aftershocks.size(); i++)
				if(aftershocks.get(i).getMag() > maxMag)
					maxMag = aftershocks.get(i).getMag();
			
			region = buildRegion(largestShock, getCentroid());
			
			// update with eventID from catalog
			mainshock.setEventId(catalogEventID);
			eventIDParam.setValue(catalogEventID);
			eventIDParam.getEditor().refreshParamEditor();
			
			// populate/validate data and forecast time windows
			updateForecastTimes();
	}

	private static double unwrap(double lon){
		if (lon > 180)
			lon -= 360;
		return lon;
	}
	
	private static String getCatalogLine(ObsEqkRupture rup) {
		StringBuilder sb = new StringBuilder();
		Location hypoLoc = rup.getHypocenterLocation();
		sb.append(catDateFormat.format(rup.getOriginTimeCal().getTime())).append("\t");
		sb.append((float)hypoLoc.getLatitude()).append("\t");
		sb.append((float)unwrap(hypoLoc.getLongitude())).append("\t");
		sb.append((float)hypoLoc.getDepth()).append("\t");
		sb.append((float)rup.getMag());
		return sb.toString();
	}
	
	private static ObsEqkRupture fromCatalogLine(String line) throws ParseException {
		line = line.trim();
		String[] split = line.split("\\s+");
		Preconditions.checkState(split.length == 10, "Unexpected number of colums. Has %s, expected 10", split.length);
		String dateStr = split[0]+"\t"+split[1]+"\t"+split[2]+"\t"+split[3]+"\t"+split[4]+"\t"+split[5];
		Date date = catDateFormat.parse(dateStr);
		double lat = Double.parseDouble(split[6]);
		double lon = Double.parseDouble(split[7]);
		double depth = Double.parseDouble(split[8]);
		double mag = Double.parseDouble(split[9]);
		Location hypoLoc = new Location(lat, lon, depth);
		
		String eventId = dateStr.replaceAll("\t", "_")+"_M"+(float)mag;
		long originTimeInMillis = date.getTime();
		
		return new ObsEqkRupture(eventId, originTimeInMillis, hypoLoc, mag);
	}
	
	private void saveCatalog(File file) throws IOException{

		FileWriter fw = null;
		try {
			workingDir = file.getParentFile();

			fw = new FileWriter(file);
			//write the header
			String headerString = "# Header: eventID " + mainshock.getEventId() + " dataEndTime " + dataEndTimeParam.getValue() + "\n";  
			fw.write(headerString);

			//write all the earthquakes
			fw.write(catalogText.getText());
			
		} catch (IOException e) {
			e.printStackTrace();
			throw e;
//			JOptionPane.showMessageDialog(this, e.getMessage(),
//					"Error Saving Catalog", JOptionPane.ERROR_MESSAGE);
		} finally {
			if (fw != null) {
				try {
					fw.close();
				} catch (IOException e) {
					e.printStackTrace();
					throw e;
				}
			}
		}
	}

	private String fromHeaderLine(String line, String target) throws ParseException {
		line = line.trim();
		String[] split = line.split("\\s+");
		
		String match = "";
		
		for (int i = 0; i < split.length; i++)
			if (split[i].equals(target)){
				if (i < split.length - 1)
					match = split[i+1];
				continue;
			}
		return match;
	}
	
	private void plotCatalogText() {
		
		StringBuilder sb = new StringBuilder();
		sb.append("# Year\tMonth\tDay\tHour\tMinute\tSec\tLat\tLon\tDepth\tMagnitude\n");
		sb.append("# Main Shock:\n");
		sb.append("# ").append(getCatalogLine(mainshock)).append("\n");
		for (ObsEqkRupture rup : aftershocks) {
			sb.append(getCatalogLine(rup)).append("\n");
		}
		if (catalogText == null) {
			catalogText = new JTextArea(sb.toString());
			catalogText.setEditable(false);
			
		} else {
			catalogText.setText(sb.toString());
		}
		
	
	}
	
	private void plotPDFs() {
		if (pdfGraphsPane == null)
			pdfGraphsPane = new JTabbedPane();
		else
			while (pdfGraphsPane.getTabCount() > 0)
				pdfGraphsPane.removeTabAt(0);
		
		HistogramFunction[] amsValExtras = null;
		HistogramFunction[] aValExtras = null;
		HistogramFunction[] pValExtras = null;
		HistogramFunction[] logcValExtras = null;
		
		if (genericModel != null) {
			if(D) System.out.println("getting generic pdfs");
			HistogramFunction genericAms = genericModel.getPriorPDF_ams(); 
			HistogramFunction genericA = genericModel.getPriorPDF_a();
			HistogramFunction genericP = genericModel.getPriorPDF_p();
			HistogramFunction genericLogC = genericModel.getPriorPDF_logc();
			
			genericAms.setName("Typical");
			genericA.setName("Typical");
			genericP.setName("Typical");
			genericLogC.setName("Typical");

			HistogramFunction bayesianAms;
			HistogramFunction bayesianA;
			HistogramFunction bayesianP;
			HistogramFunction bayesianLogC;
			if (bayesianModel != null) {
				
				if(D) System.out.println("getting bayesian pdfs");
				bayesianAms = bayesianModel.getPDF_ams(); 
				bayesianA = bayesianModel.getPDF_a();
				bayesianP = bayesianModel.getPDF_p();
				bayesianLogC = bayesianModel.getPDF_logc();
				
				if(bayesianAms != null)
					bayesianAms.setName("This sequence");
				if(bayesianA != null)
					bayesianA.setName("This sequence");
				if(bayesianP != null)
					bayesianP.setName("This sequence");
				if(bayesianLogC != null)
					bayesianLogC.setName("This sequence");
				
				amsValExtras = new HistogramFunction[] { genericAms, bayesianAms };
				aValExtras = new HistogramFunction[] { genericA, bayesianA };
				pValExtras = new HistogramFunction[] { genericP, bayesianP };
				logcValExtras = new HistogramFunction[] { genericLogC, bayesianLogC };
			} else {
				amsValExtras = new HistogramFunction[] { genericAms };
				aValExtras = new HistogramFunction[] { genericA };
				pValExtras = new HistogramFunction[] { genericP };
				logcValExtras = new HistogramFunction[] { genericLogC };
			}
		}
		
		if(D) System.out.println("getting sequence specific pdfs");
		add1D_PDF(bayesianModel.getPDF_ams(), "ams-value", amsValExtras);	
		add1D_PDF(bayesianModel.getPDF_a(), "a-value", aValExtras);
		add1D_PDF(bayesianModel.getPDF_p(), "p-value", pValExtras);
		add1D_PDF(bayesianModel.getPDF_logc(), "logc-value", logcValExtras);
		add2D_PDF(bayesianModel.get2D_PDF_for_a_and_logc(), "a-value", "logc-value");
		add2D_PDF(bayesianModel.get2D_PDF_for_a_and_p(), "a-value", "p-value");
		add2D_PDF(bayesianModel.get2D_PDF_for_logc_and_p(), "logc-value", "p-value");
		
		Runnable displayRun = new Runnable() {
			
			@Override
			public void run() {
				if (tabbedPane.getTabCount() >= aftershock_expected_index){ //&& !tabbedPane.isEnabledAt(aftershock_expected_index)){
					tabbedPane.removeTabAt(pdf_tab_index);
					tabbedPane.insertTab("Model PDFs", null, pdfGraphsPane,
							"Aftershock Model Prob Dist Funcs", pdf_tab_index);
				}
				else
					Preconditions.checkState(tabbedPane.getTabCount() > pdf_tab_index, "Plots added out of order");

			}
		};
		
		if (SwingUtilities.isEventDispatchThread()) {
			displayRun.run();
		} else {
			try {
				SwingUtilities.invokeAndWait(displayRun);
			} catch (Exception e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
	}
	
	private void add1D_PDF(HistogramFunction pdf, String name, HistogramFunction... extras) {
		if (pdf == null)
			return;
		
		Preconditions.checkState(!Double.isNaN(pdf.getMaxY()), "NaN found in "+pdf.getName());
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(pdf);
		List<PlotCurveCharacterstics> chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, sequence_specific_color));
		
		if (extras != null && extras.length > 0) {
			for (int i=0; i<extras.length; i++) {
				funcs.add(extras[i]);
				Color c;
				String extraName = extras[i].getName().toLowerCase();
				if (extraName.contains("generic"))
					c = generic_color;
				else if (extraName.contains("bayesian"))
					c = bayesian_color;
				else
					c = extra_colors[i % extra_colors.length];
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
			}
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, pdf.getName(), name, "Density");
		spec.setLegendVisible(funcs.size() > 1);

		GraphWidget widget = new GraphWidget(spec);
		setupGP(widget);
		widget.getButtonControlPanel().setVisible(false);
		GraphPanel graphPanel = widget.getGraphPanel();
		graphPanel.getComponent(2).setVisible(false);
		
		pdfGraphsPane.addTab(name, null, widget);
		
		
	}
	
	private void add2D_PDF(EvenlyDiscrXYZ_DataSet pdf, String name1, String name2) {
		if (pdf == null)
			return;
		
		String title = "PDF for "+name1+" vs "+name2;
		
		Preconditions.checkState(Doubles.isFinite(pdf.getMaxZ()), "NaN found in "+title);
		
		CPT cpt = getDistCPT().rescale(pdf.getMinZ(), pdf.getMaxZ());
		
		XYZPlotSpec spec = new XYZPlotSpec(pdf, cpt, title, name1, name2, "Density");
		PlotPreferences prefs = XYZGraphPanel.getDefaultPrefs();
	
		XYZGraphPanel xyzGP = new XYZGraphPanel(prefs);
		xyzGP.setPreferredSize(new Dimension(chartWidth, chartHeight));
		
		pdfGraphsPane.addTab(name1+" vs "+name2, null, xyzGP);
		double xDelta = pdf.getGridSpacingX();
		double yDelta = pdf.getGridSpacingY();
		xyzGP.drawPlot(spec, false, false,
				new Range(pdf.getMinX()-0.5*xDelta, pdf.getMaxX()+0.5*xDelta),
				new Range(pdf.getMinY()-0.5*yDelta, pdf.getMaxY()+0.5*yDelta));
		
		Component graphComponent = xyzGP.getComponent(0);
		graphComponent.setPreferredSize(new Dimension(chartWidth-100,chartHeight-60-100));
		
	}
	
	private class CalcStep {
		
		private String title;
		private String progressMessage;
		private Runnable run;
		private boolean runInEDT;
		
		public CalcStep(String title, String progressMessage, Runnable run) {
			this(title, progressMessage, run, false);
		}
		
		public CalcStep(String title, String progressMessage, Runnable run, boolean runInEDT) {
			this.title = title;
			this.progressMessage = progressMessage;
			this.run = run;
			this.runInEDT = runInEDT;
			
		}
	}
	
	private volatile boolean stopRequested = false;	//this should all be revised to use an Executor or something
	private volatile boolean pauseRequested = false;	//but this works for now at keeping UX pretty smooth
	private class CalcRunnable implements Runnable {
		private CalcProgressBar progress;
		private CalcStep[] steps;
		private String curTitle;
		private Throwable exception;
		
		public CalcRunnable(CalcProgressBar progress, CalcStep... calcSteps) {
			this.progress = progress;
			this.steps = calcSteps;
		}
		
		@Override
		public void run() {
			
			pauseRequested = true;
			SwingUtilities.invokeLater(new Runnable() {
			
				@Override
				public void run() {
					if (progress.isVisible())
						progress.setVisible(false);
					
					progress.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
					progress.setModalityType(ModalityType.APPLICATION_MODAL);
					
					
					if(steps.length > 1){
						progress.updateProgress(0, steps.length);
					}else
						progress.setIndeterminate(true);
					
					try{
						progress.pack();
						progress.setVisible(true);
												
					} catch(Exception e) {
						System.err.println("Afershock forecaster has encountered an error and needs to close.\n"
								+ "Please restart the aftershock forecasting software. Unpredictable behavior may occur without restart.\n");
					} finally {
						pauseRequested = false;
					}
				}
			});
			
			WindowListener wl = new WindowAdapter() {
				@Override
				public void windowClosing(WindowEvent e){
					if(verbose){
						System.out.println("Calculation interrupted at step: " + steps[0].progressMessage);
						System.out.println("remaining steps:");
						for (int i = 1; i < steps.length; i++){
							if (0 < steps.length)
								System.out.println("\t"+ steps[i].progressMessage);
						}
					}
					stopRequested = true;
				}
			};
			progress.addWindowListener(wl);
			
			int stepCount = -1;
			stopRequested = false;
			setChangeListenerEnabled(false); //no cascading updates. All routines must be initiated from EDT, not the changeListener
			while(!stopRequested && ++stepCount < steps.length){
				
				progress.removeWindowListener(wl);
				final int sc = stepCount; 
				wl = new WindowAdapter() {
					
					@Override
					public void windowClosing(WindowEvent e){
						if(verbose){
							System.out.println("Calculation interrupted at step: " + steps[sc].progressMessage);
							System.out.println("remaining steps:");
							for (int i = sc+1; i < steps.length; i++){
								if (sc < steps.length)
									System.out.println("\t"+ steps[i].progressMessage);
							}
						}
						stopRequested = true;
					}
				};
				progress.addWindowListener(wl);
				
				final CalcStep step = steps[stepCount];
				final int innerStepCount = stepCount; 

				pauseRequested = true;
				SwingUtilities.invokeLater(new Runnable() {

					@Override
					public void run() {
						progress.setTitle(step.title);
						progress.updateProgress(innerStepCount+1, steps.length, step.progressMessage);
						progress.setProgressMessage(step.progressMessage);
						progress.pack();
						pauseRequested = false;
					}
				});

				long timeout = 1000;
				long timesofar = 0;
				while (pauseRequested && timesofar < timeout) {
					try {
						if(D) System.out.println("...");
						timesofar += 100;
						Thread.sleep(100);
						
					} catch (InterruptedException e1) {
						// TODO Auto-generated catch block
						System.err.println("Calculation \"" + step.title + "\" was interrupted");
					}
				}
				
				
				curTitle = step.title;
				if (step.runInEDT) {
					
					try {
						SwingUtilities.invokeAndWait(new Runnable() {

							@Override
							public void run() {
								
									try {
										step.run.run();
										pauseRequested = false;
									} catch (Throwable e) {
										exception = e;
									}
								
							}
						});
					} catch (Exception e) {
						exception = e;
					}
					if (exception != null)
						break;
					
//					pauseRequested = true;
//					SwingUtilities.invokeLater(new Runnable() {
//
//						@Override
//						public void run() {
//							try {
//								step.run.run();
//								pauseRequested = false;
//							} catch (Throwable e) {
//								exception = e;
//							}
//						}
//					});

				} else {
					try {
						step.run.run();
						pauseRequested = false;
					} catch (Throwable e) {
						exception = e;
						break;
					}
				}
			}

			while (pauseRequested){
				if(D) System.out.println("...");
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}
			
			
			setChangeListenerEnabled(true); //ok, done. start listening again.

				SwingUtilities.invokeLater(new Runnable() {
					@Override
					public void run() {
						progress.setVisible(false);
						progress.dispose();
					}
				});
			
			if (exception != null) {
				final String title = "Error "+curTitle;
				exception.printStackTrace();
				final String message = exception.getMessage();
				try {
					SwingUtilities.invokeAndWait(new Runnable() {
						
						@Override
						public void run() {
							JOptionPane.showMessageDialog(null, message, title, JOptionPane.ERROR_MESSAGE);
						}
					});
				} catch (Exception e) {
					System.err.println("Error displaying error message!");
					e.printStackTrace();
				}
				
			}
		}
	}

	@Override
	public void parameterChange(ParameterChangeEvent event) {
		Parameter<?> param = event.getParameter();
		
		final CalcProgressBar progress = new CalcProgressBar(this, "Progress", "Intializing...", false);
		
		// define common steps
		CalcStep fetchCatalogStep = new CalcStep("Fetching Events",
				"Contacting USGS ComCat Webservice. This is occasionally slow. "
						+ "If it fails, trying again often works.", new Runnable() {

			@Override
			public void run() {
				// set the data window to 0 - now, forecast time to now - one month
				dataStartTimeParam.setValue(0);
				if (dataEndTimeParam.getValue() == null)
					dataEndTimeParam.setValue(dataEndTimeParam.getMax());
				
				fetchEvents();
				if (mainshock == null) stopRequested=true;
			}
			
		}, true);

		CalcStep setMainshockStep = new CalcStep("Initializing Model", "Collecting mainshock information...", new Runnable() {

			@Override
			public void run() {
				setMainshock(mainshock);
			}
			
		}, true);

		CalcStep postFetchCalcStep = new CalcStep("Initializing Model", "Computing sequence statistics...", new Runnable() {

			@Override
			public void run() {
				doPostFetchCalculations();
			}
		}, true);
		
//		computeBButton
		CalcStep computeBStep = new CalcStep("Initializing Model", "Computing b-value and magnitude of completeness...", new Runnable() {
			
			@Override
			public void run() {

				if (aftershocks.size() > 0){
					ogataMND = new OgataMagFreqDist(aftershocks, genericModel.get_b(), genericModel.get_bSigma());

					if(verbose) System.out.println(ogataMND + " " + ogataMND.goodnessOfFit());
				
					Double b, mc;
					b = ogataMND.get_bValue();
					mc = ogataMND.get_Mc();
					validateParameter(mc, "Mc");
				
					bParam.setValue(round(b,2));
					bParam.getEditor().refreshParamEditor();

					mcParam.setValue(round(mc,2));
					mcParam.getEditor().refreshParamEditor();

					autoMcParam.setValue(false);
					computeBButton.setButtonText("Reset Mc");
					computeBButton.setInfo("Estimate Mc from b-value stability");
					computeBButton.getEditor().refreshParamEditor();
					
//					if(seqSpecModel != null)
//						seqSpecModel.set_b(b);
					if(bayesianModel != null)
						bayesianModel.set_b(b);
					
					link_alpha(); //reset the alpha parameter, too
					
//					tabbedPane.setSelectedIndex(mag_num_tab_index);
				} else {
					System.err.println("Not enough data to compute b-value");
				}
				
				setEnableParamsPostComputeB(true);
				
			}
			
		}, true);

		CalcStep postFetchPlotStep = new CalcStep("Plotting", "Building sequence summary plots...", new Runnable() {

			@Override
			public void run() {
				doPostFetchPlots();
				//inlcudes MFD plot
			}
		}, true);
		
		CalcStep updateModelStep = new CalcStep("Updating Model", "Updating generic model parameters...", new Runnable() {

			@Override
			public void run() {
				// if a generic model has already been computed, update it to reflect new Mc/Mref productivity adjustment
				if(genericModel != null){
					if(verbose) System.out.println("Adjusting generic productivity to reflect Mc change...");

					// not sure if I need this
					genericModel = new ETAS_AftershockModel_Generic(mainshock, aftershocks, genericParams, 
							dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
							forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), mcParam.getValue(),
							genericParams.get_maxMag(), 100, 0, fitMSProductivityParam.getValue(), timeDepMcParam.getValue());

					resetFitConstraints(genericParams);
					setEnableParamsPostAftershockParams(false);
				}
			}


		}, true);
		
//		computeAftershockParamsButton
//		CalcStep computeAftershockParamStep = new CalcStep("Estimating ETAS parameters", "Computing sequence-specific model. This may take some time...", new Runnable() {
//
//			@Override
//			public void run() {
//				Range amsRange = amsValRangeParam.getValue();
//				int amsNum = amsValNumParam.getValue();
//				validateRange(amsRange, amsNum, "ams-value");
//				
//				Range aRange = aValRangeParam.getValue();
//				int aNum = aValNumParam.getValue();
//				validateRange(aRange, aNum, "a-value");
//				
//				Range pRange = pValRangeParam.getValue();
//				int pNum = pValNumParam.getValue();
//				validateRange(pRange, pNum, "p-value");
//				
//				Range cRange = cValRangeParam.getValue();
//				int cNum = cValNumParam.getValue();
//				validateRange(cRange, cNum, "c-value");
//				
//				Double mc = mcParam.getValue();
//				validateParameter(mc, "Mc");
//								
//				Double b = bParam.getValue();
//				validateParameter(b, "b-value");
//				
//				Double alpha = alphaParam.getValue();
//								
//				if (timeDepMcParam.getValue()) {
//					
//					//to do: replace with rmax routine from desktop version of the code
//					seqSpecModel = new ETAS_AftershockModel_SequenceSpecific(mainshock, aftershocks,
//							ETAS_StatsCalc.linspace(amsRange.getLowerBound(), amsRange.getUpperBound(), amsNum), 
//							genericParams.get_aSigma(),
//							ETAS_StatsCalc.linspace(aRange.getLowerBound(), aRange.getUpperBound(), aNum),
//							ETAS_StatsCalc.linspace(pRange.getLowerBound(), pRange.getUpperBound(), pNum),
//							ETAS_StatsCalc.logspace(cRange.getLowerBound(), cRange.getUpperBound(), cNum),
//							alpha, b, magRefParam.getValue(),	dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
//							forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(),	// forecast start time, forecast end time 
//							mcParam.getValue(), genericParams.get_maxMag(), 100, 0, fitMSProductivityParam.getValue(), timeDepMcParam.getValue(),
////							new ETAS_AftershockModel_Generic(genericParams), progress);	//max sim mag, max number of generations, number of simulations
//							genericModel, progress);	//max sim mag, max number of generations, number of simulations
//					
//					
//				} else {
//					seqSpecModel = new ETAS_AftershockModel_SequenceSpecific(mainshock, aftershocks, 
//							ETAS_StatsCalc.linspace(amsRange.getLowerBound(), amsRange.getUpperBound(), amsNum), 
//							genericParams.get_aSigma(),
//							ETAS_StatsCalc.linspace(aRange.getLowerBound(), aRange.getUpperBound(), aNum),
//							ETAS_StatsCalc.linspace(pRange.getLowerBound(), pRange.getUpperBound(), pNum),
//							ETAS_StatsCalc.logspace(cRange.getLowerBound(), cRange.getUpperBound(), cNum),
//							alpha, b, magRefParam.getValue(), dataStartTimeParam.getValue(), dataEndTimeParam.getValue(), 
//							forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(),
//							mcParam.getValue(), genericParams.get_maxMag(), 100, 0, fitMSProductivityParam.getValue(), timeDepMcParam.getValue(),
//							null, progress);
//							
//				}
//				if(verbose) System.out.format("Mainshock productivity magnitude: %2.2f\n" , seqSpecModel.getProductivityMag());
//			}
//			
//		}, true);
		
		CalcStep computeBayesStep = new CalcStep("Estimating ETAS parameters", "Computing Bayesian model. This may take some time...", new Runnable() {

			@Override
			public void run() {

//				bayesianModel = null;
//				if (genericModel != null) {
//					bayesianModel = new ETAS_AftershockModel_SequenceSpecific(mainshock, aftershocks, rmaxParam.getValue(),
//							seqSpecModel.ams_vec, genericParams.get_aSigma(), seqSpecModel.a_vec, seqSpecModel.p_vec,	seqSpecModel.c_vec,
//							seqSpecModel.alpha, seqSpecModel.b, magRefParam.getValue(), dataStartTimeParam.getValue(), dataEndTimeParam.getValue(), 
//							forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(),
//							mcParam.getValue(), genericParams.get_maxMag(), 100, 0, fitMSProductivityParam.getValue(), timeDepMcParam.getValue(),
//							genericModel, progress);
//				}
				Range amsRange = amsValRangeParam.getValue();
				int amsNum = amsValNumParam.getValue();
				validateRange(amsRange, amsNum, "ams-value");
				
				Range aRange = aValRangeParam.getValue();
				int aNum = aValNumParam.getValue();
				validateRange(aRange, aNum, "a-value");
				
				Range pRange = pValRangeParam.getValue();
				int pNum = pValNumParam.getValue();
				validateRange(pRange, pNum, "p-value");
				
				Range cRange = cValRangeParam.getValue();
				int cNum = cValNumParam.getValue();
				validateRange(cRange, cNum, "c-value");
				
				Double mc = mcParam.getValue();
				validateParameter(mc, "Mc");
								
				Double b = bParam.getValue();
				validateParameter(b, "b-value");
				
				Double alpha = alphaParam.getValue();
				bayesianModel = null;
				if (genericModel != null) {
					bayesianModel = new ETAS_AftershockModel_SequenceSpecific(mainshock, aftershocks, 
							ETAS_StatsCalc.linspace(amsRange.getLowerBound(), amsRange.getUpperBound(), amsNum), 
							genericParams.get_aSigma(),
							ETAS_StatsCalc.linspace(aRange.getLowerBound(), aRange.getUpperBound(), aNum),
							ETAS_StatsCalc.linspace(pRange.getLowerBound(), pRange.getUpperBound(), pNum),
							ETAS_StatsCalc.logspace(cRange.getLowerBound(), cRange.getUpperBound(), cNum),
							alpha, b, magRefParam.getValue(),	dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
							forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(),
							mcParam.getValue(), genericParams.get_maxMag(), 100, 0, fitMSProductivityParam.getValue(), timeDepMcParam.getValue(),
							genericModel, progress);
				}


				amsValParam.setValue(round(bayesianModel.getMaxLikelihood_ams(),2));
				amsValParam.getEditor().refreshParamEditor();
				aValParam.setValue(round(bayesianModel.getMaxLikelihood_a(),2));
				aValParam.getEditor().refreshParamEditor();
				pValParam.setValue(round(bayesianModel.getMaxLikelihood_p(),2));
				pValParam.getEditor().refreshParamEditor();
//				cValParam.setValue(round(bayesianModel.getMaxLikelihood_c(),6));
				cValParam.setValue(round(Math.log10(bayesianModel.getMaxLikelihood_c()),2));
				cValParam.getEditor().refreshParamEditor();
			}
		}, true);
		
		CalcStep postParamPlotStep = new CalcStep("Plotting", "Building ETAS summary plots...", new Runnable() {
			
			@Override
			public void run() {
	
				plotPDFs();
				setEnableParamsPostAftershockParams(true);
				plotCumulativeNum();

			}
		}, true);
		
//		computeAftershockForecastButton
		CalcStep genericForecastStep = new CalcStep("Computing aftershock forecast", "Generating Generic forecast. This can take some time...",
				new Runnable() {

					@Override
					public void run() {
						if(verbose) System.out.println("Computing Generic forecast");
						genericModel.setMagComplete(mcParam.getValue());
						genericModel.computeNewForecast(dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
								forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), 10000);
					}
		}, true);
		
//		CalcStep seqSpecForecastStep = new CalcStep("Computing aftershock forecast", "Generating Sequence-Specific forecast. This can take some time...",
//				new Runnable() {
//
//					@Override
//					public void run() {
//
//						if (seqSpecModel != null){
//							seqSpecModel.setMagComplete(mcParam.getValue());
//							if(verbose)System.out.println("Computing sequence-only forecast");
//							seqSpecModel.computeNewForecast(dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
//									forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), 10000);
//						}
//					}
//		}, true);
		
		CalcStep bayesianForecastStep = new CalcStep("Computing aftershock forecast", "Generating Bayesian forecast. This can take some time...",
				new Runnable() {

					@Override
					public void run() {
						
						if (bayesianModel != null){
							bayesianModel.setMagComplete(mcParam.getValue());
							if(verbose)System.out.println("Computing Bayesian forecast");
							bayesianModel.computeNewForecast(dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
									forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), 10000);
						}
					}		
		}, true);
		
		CalcStep plotMFDStep = new CalcStep("Generating Forecast Summary", "Building summary plots...", new Runnable() {

			@Override
			public void run() {
				if (tabbedPane.getTabCount() > mag_time_tab_index)
					plotMFDs(aftershockMND, mcParam.getValue());
			}
		}, true);
		
		CalcStep postForecastPlotStep = new CalcStep("Generating Forecast Summary", "Building forecast plots...",
				new Runnable() {
					@Override
					public void run() {
						plotExpectedAfershockMFDs();
						//update cumNumPlot
						plotCumulativeNum();
						setEnableParamsPostForecast(true);
					}
			
		}, true);
		
		CalcStep plotTableStep = new CalcStep("Generating Forecast Summary", "Building aftershock forecast table...",
				new Runnable() {

					@Override
					public void run() {
						plotForecastTable(progress);
//						tabbedPane.setSelectedIndex(aftershock_expected_index);
					}
		}, true);	
		
//		generateMapButton
		CalcStep fetchSourceStep = new CalcStep("Fetching Finite Source", "Contacting USGS ComCat Webservice. If this doesn't work, trying again may help...",
				new Runnable() {
					@Override
					
					public void run(){
						
						if(verbose) System.out.println("Fetching ShakeMap finite source from ComCat...");
//						if (!fitMainshockSourceParam.getValue() && faultTrace == null){
						final FitSourceType fitType = fitSourceTypeParam.getValue();
						
						if (fitType == FitSourceType.SHAKEMAP){
								faultTrace = accessor.fetchShakemapSource(eventIDParam.getValue());
								if(faultTrace == null){
									System.out.println("...ShakeMap finite source not available. Fitting early aftershocks...");
									// no shakemap, fit early aftershocks
									fitSourceTypeParam.setValue(FitSourceType.AFTERSHOCKS);
									fitSourceTypeParam.getEditor().refreshParamEditor();
									faultTrace = new FaultTrace("Mainshock Location");
									faultTrace.add(mainshock.getHypocenterLocation());
								}
								
								if (mainshock.getHypocenterLocation().getLongitude() > 180){
									ListIterator<Location> i = faultTrace.listIterator();
									while (i.hasNext()){
										Location loc = i.next();
										if (loc.getLongitude() < 0){
											i.set(new Location(loc.getLatitude(), loc.getLongitude()+360));
										}
									}
								}
								if(D) System.out.println(faultTrace);
						}
						
						if(verbose) System.out.println("Fetching ShakeMap URL from ComCat servers...");
						shakeMapURL = accessor.fetchShakemapURL(eventIDParam.getValue());
						if(verbose) System.out.println(shakeMapURL);
							
						
					}
		}, true);
		
		CalcStep plotMapStep = new CalcStep("Generating Spatial Forecast", "Generating spatial forecast maps. This can take some time...",
				new Runnable() {

					@Override
					public void run() {
						Stopwatch watch = Stopwatch.createStarted();

						plotRateModel2D(progress);
						setEnableParamsPostRender(true);
						
						watch.stop();
						if(verbose)System.out.println(
								"Took "+watch.elapsed(TimeUnit.SECONDS)+"s to compute/plot rate map");
						//								tabbedPane.setSelectedIndex(aftershock_expected_index);
						tabbedPane.setSelectedIndex(forecast_map_tab_index);
						
						
					}
					
		}, true);
		
		CalcStep tabRefreshStep = new CalcStep("Plotting","Updating plots...", new Runnable() {
			
			@Override
			public void run() {
				refreshTabs(forecast_map_tab_index);
			}
			
		}, true);
		
		CalcStep tabSelectStep = new CalcStep("Plotting","Updating plots...", new Runnable() {
			
			@Override
			public void run() {
				tabbedPane.setSelectedIndex(forecast_table_tab_index);
			}
			
		}, true);
		
		CalcStep autoSaveStep = new CalcStep("Saving catalog","Saving catalog to local file...",
				new Runnable() {
					@Override
					public void run() {
						String filename = new String(workingDir + "/Aftershock_Catalog.txt");
						File savefile = new File(filename);
						try {
							saveCatalog(savefile);
						} catch (IOException e) {
							ExceptionUtils.throwAsRuntimeException(e);
						}
					}
				}, true);
		

		CalcStep tipStep; //define later

		if (!changeListenerEnabled || param.getValue() == null) {
			if(D)	System.out.println("Suppressing refresh for " + param.getName());
		} else {
			if(verbose) System.out.println("Updating " + param.getName());
			// putting this in as a safegaurd
			setChangeListenerEnabled(false);
			try {
			if (param == quickForecastButton) {
				System.out.println("Computing quick forecast with defafult settings...");

				CalcStep quickTabRefreshStep1 = new CalcStep("Plotting","Updating plots...", new Runnable() {
					@Override
					public void run() {
						refreshTabs(cml_num_tab_index);
					}
				}, true);
				CalcStep quickTabRefreshStep2 = new CalcStep("Plotting","Updating plots...", new Runnable() {
					@Override
					public void run() {
						refreshTabs(cml_num_tab_index);
					}
				}, true);
				CalcStep quickTabRefreshStep3 = new CalcStep("Plotting","Updating plots...", new Runnable() {
					@Override
					public void run() {
						refreshTabs(aftershock_expected_index);
					}
				}, true);
				CalcStep quickTabRefreshStep4 = new CalcStep("Plotting","Updating plots...", new Runnable() {
					@Override
					public void run() {
						refreshTabs(forecast_map_tab_index);
					}
				}, true);


				CalcStep quickTabSelectStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
						tabbedPane.setSelectedIndex(forecast_map_tab_index);
					}

				}, true);

				tipStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
						printTip(6);
					}

				}, true);

				CalcRunnable run = new CalcRunnable(progress,
						fetchCatalogStep,
						setMainshockStep,
						postFetchCalcStep,
						postFetchPlotStep,
//						computeBStep, //not stable enough...
						plotMFDStep,
						quickTabRefreshStep1,
//						computeAftershockParamStep,
						computeBayesStep,
						postParamPlotStep,
						quickTabRefreshStep2,
						genericForecastStep,
//						seqSpecForecastStep,
						bayesianForecastStep,
						postForecastPlotStep,
						plotTableStep,
						quickTabRefreshStep3,
						fetchSourceStep,
						plotMapStep,
						quickTabRefreshStep4,
						quickTabSelectStep,
						tipStep);
				new Thread(run).start();

			} else if (param == eventIDParam || param == dataStartTimeParam || param == dataEndTimeParam
					|| param == regionEditParam) {
//				if(verbose) System.out.println("Updating " + param.getName()); 

				if (param == regionEditParam){
					// check for international dateline crossing
					try{
						if (minLonParam.getValue() > 0 && maxLonParam.getValue() < 0)
							wrapLongitudes(true);
					} catch (Exception e){
						//doesn't exist
					}
				}
				
				// reset the workspace
				resetWorkspace();

				if(param == eventIDParam) printTip(1);

			} else if (param == nowBoolean){
//				if(verbose) System.out.println("Updating start now flag");
				refreshTimeWindowEditor();

			} else if (param == forecastDurationParam) {
//				if(verbose) System.out.println("Updating forecast duration");
				try{
					if (forecastStartTimeParam.getValue() != null){
						forecastEndTimeParam.setValue(round(forecastStartTimeParam.getValue() + forecastDurationParam.getValue()));
						forecastEndTimeParam.getEditor().refreshParamEditor();
					}
				}finally{

				}
			
			} else if (param == forecastStartTimeParam) {
//				if(verbose) System.out.println("Updating forecast start time");
				setEnableParamsPostFetch(false);

				try{
					dataEndTimeParam.setValue(forecastStartTimeParam.getValue());
					dataEndTimeParam.getEditor().refreshParamEditor();	
					forecastEndTimeParam.setValue(forecastStartTimeParam.getValue() + forecastDurationParam.getValue());
					forecastEndTimeParam.getEditor().refreshParamEditor();
				} finally {
				}
				refreshTimeWindowEditor();

			} else if (param == regionTypeParam || param == regionCenterTypeParam) {
//				if(verbose) System.out.println("Updating regionType");
				updateRegionParamList(regionTypeParam.getValue(), regionCenterTypeParam.getValue());
		
			} else if (param == fetchButton) {
				System.out.println("Fetching events from ComCat...");
				setEnableParamsPostFetch(false);

				CalcStep fetchTabRefreshStep = new CalcStep("Plotting",  "Updating plots...", new Runnable() {

					@Override
					public void run() {
						refreshTabs(cml_num_tab_index);
					}
				}, true);

				CalcStep fetchTabSelectStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
						tabbedPane.setSelectedIndex(epicenter_tab_index);
					}

				}, true);

				tipStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
						printTip(2);
					}

				}, true);

				CalcRunnable run = new CalcRunnable(progress,
						fetchCatalogStep, 
						setMainshockStep,
						postFetchCalcStep,
						postFetchPlotStep,
						fetchTabRefreshStep,
						fetchTabSelectStep,
						tipStep);
				new Thread(run).start();

			} else if (param == loadCatalogButton) {
				System.out.println("Loading Catalog from file...");

				if (loadCatalogChooser == null)
					loadCatalogChooser = new JFileChooser();

				loadCatalogChooser.setCurrentDirectory(workingDir);

				int ret = loadCatalogChooser.showOpenDialog(this);
				if (ret == JFileChooser.APPROVE_OPTION) {
					final File file = loadCatalogChooser.getSelectedFile();

					CalcStep loadStep = new CalcStep("Loading Catalog", "Loading events from local catalog...", new Runnable() {

						@Override
						public void run() {
							try {
								loadCatalog(file);
								setEnableParamsPostFetch(false);

							} catch (IOException e) {
								ExceptionUtils.throwAsRuntimeException(e);
							}
						}
					}, true);


					tabRefreshStep = new CalcStep("Plotting",  "Updating plots...", new Runnable() {

						@Override
						public void run() {
							refreshTabs(cml_num_tab_index);
						}
					}, true);

					tabSelectStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

						@Override
						public void run() {
							tabbedPane.setSelectedIndex(epicenter_tab_index);
						}

					}, true);


					tipStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

						@Override
						public void run() {
							printTip(2);
						}

					}, true);

					CalcRunnable run = new CalcRunnable(progress, loadStep, setMainshockStep, postFetchCalcStep, postFetchPlotStep,
							tabRefreshStep, tabSelectStep, tipStep);
					new Thread(run).start();

				}

			} else if (param == saveCatalogButton) {
				System.out.println("Saving Catalog to file...");
				
				if (saveCatalogChooser == null)
					saveCatalogChooser = new JFileChooser();

				saveCatalogChooser.setCurrentDirectory(workingDir);

				int ret = saveCatalogChooser.showSaveDialog(this);
				if (ret == JFileChooser.APPROVE_OPTION) {
					File file = saveCatalogChooser.getSelectedFile();
					
					CalcStep saveStep = new CalcStep("Saving Catalog", "Saving catalog to file...", new Runnable() {

						@Override
						public void run() {
							try {
								saveCatalog(file);
							} catch (IOException e) {
								ExceptionUtils.throwAsRuntimeException(e);
							}
						}
					}, true);
					CalcRunnable run = new CalcRunnable(progress, saveStep);
					new Thread(run).start();
					
				}
				
			} else if (param == mcParam){ // || param == magPrecisionParam) {
//				if(verbose) System.out.println("Updating Mc");
				
				changeListenerEnabled = false;
				
				setEnableParamsPostComputeMc(false);

				// check that it's not greater than mainshock magnitude
				if(mcParam.getValue() > mainshock.getMag() - 0.05){
					System.err.println("Mc cannot be set to larger than the mainshock magnitude");
					mcParam.setValue(round(mainshock.getMag() - 0.05,2));
					mcParam.getEditor().refreshParamEditor();
				} 
				
				changeListenerEnabled = true;

				CalcStep setMcStep = new CalcStep("Calculating","Recalculating Mc with generic b-value...", new Runnable() {

					@Override
					public void run() {
						setMagComplete(mcParam.getValue());
						
						autoMcParam.setValue(false);
						computeBButton.setButtonText("Reset Mc");
						computeBButton.setInfo("Estimate Mc from b-value stability");
						computeBButton.getEditor().refreshParamEditor();
					}
					
				}, true);
				
				CalcStep plotStep = new CalcStep("Plotting", "Plotting magnitude-frequency distribution...", new Runnable() {

					@Override
					public void run() {
						if (tabbedPane.getTabCount() > mag_time_tab_index)
							plotMFDs(aftershockMND, mcParam.getValue());
						
						if (tabbedPane.getTabCount() > cml_num_tab_index)
							plotCumulativeNum();
					}
				}, true);

				CalcRunnable run = new CalcRunnable(progress, setMcStep, updateModelStep, plotStep);
				new Thread(run).start();

			} else if (param == computeBButton) {
				System.out.println("Calculating b-value and Mc...");

				CalcStep setMcStep = new CalcStep("Calculating","Recalculating Mc with generic b-value...", new Runnable() {

					@Override
					public void run() {
						setMagComplete();
					}
					
				}, true);
				
				tabRefreshStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
						if (tabbedPane.getTabCount() > cml_num_tab_index)
							plotCumulativeNum();
						refreshTabs(mag_num_tab_index);
					}

				}, true);

				tabSelectStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
						tabbedPane.setSelectedIndex(mag_num_tab_index);
					}

				}, true);

				CalcRunnable run;
				if (autoMcParam.getValue()){
					run = new CalcRunnable(progress, computeBStep, updateModelStep, plotMFDStep, tabRefreshStep, tabSelectStep);
					autoMcParam.setValue(false);
					computeBButton.setButtonText("Reset Mc");
					computeBButton.setInfo("Estimate Mc from b-value stability");
					computeBButton.getEditor().refreshParamEditor();
				}else{
					run = new CalcRunnable(progress, setMcStep, updateModelStep, plotMFDStep, tabRefreshStep, tabSelectStep);
					autoMcParam.setValue(true);
					computeBButton.setButtonText("Compute b");
					computeBButton.setInfo("Simultaneously estimate b and Mc.");
					computeBButton.getEditor().refreshParamEditor();
				}
				
				new Thread(run).start();
	
			} else if (param == bParam) {
				CalcStep linkAlphaStep = new CalcStep("Updating", "Updating Parameters...", new Runnable() {
					@Override
					public void run() {
						link_alpha();
					}
				}, true);

				CalcStep plotStep = new CalcStep("Plotting", "Building magnitude-frequency distribution plot...", new Runnable() {

					@Override
					public void run() {
						if (tabbedPane.getTabCount() > mag_time_tab_index)
							plotMFDs(aftershockMND, mcParam.getValue());
					}
				}, true);

				CalcRunnable run = new CalcRunnable(progress, linkAlphaStep, plotStep);
				new Thread(run).start();

			} else if (param == fitMSProductivityParam) {

				CalcStep updateParamsStep = new CalcStep("Updating", "Updating Parameters...", new Runnable() {
					@Override
					public void run() {
						//grey out the values
						if(fitMSProductivityParam.getValue() != null){
							amsValRangeParam.getEditor().setEnabled(fitMSProductivityParam.getValue());
							amsValNumParam.getEditor().setEnabled(fitMSProductivityParam.getValue());
						}
						setEnableParamsPostAftershockParams(false);
					}
				}, true);
				
				CalcRunnable run = new CalcRunnable(progress, updateParamsStep);
				new Thread(run).start();
				
			} else if (param == tectonicRegimeParam) {

				CalcStep updateTectonicRegimeStep = new CalcStep("Updating Model", "Updating generic parameters to " + tectonicRegimeParam.getValue() + "...", new Runnable() {

					@Override
					public void run() {
						setEnableParamsPostAftershockParams(false);
						setEnableParamsPostComputeB(false);

						genericParams = genericFetch.get(tectonicRegimeParam.getValue());
						Preconditions.checkNotNull(genericParams, "Generic params not found or server error");
						if(verbose) System.out.println("Updating Generic model params to " + tectonicRegimeParam.getValue() + ": "+genericParams);

						setMagComplete();
					}
				}, true);

				CalcRunnable run = new CalcRunnable(progress, updateTectonicRegimeStep, updateModelStep);
				new Thread(run).start();

			} else if (param == amsValRangeParam || param == amsValNumParam) {
				
				updateRangeParams(amsValRangeParam, amsValNumParam, 51);
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == aValRangeParam || param == aValNumParam) {
				
				updateRangeParams(aValRangeParam, aValNumParam, 51);
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == pValRangeParam || param == pValNumParam) {
				
				updateRangeParams(pValRangeParam, pValNumParam, 51);
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == cValRangeParam || param == cValNumParam) {
				
				updateRangeParams(cValRangeParam, cValNumParam, 51);
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == timeDepMcParam) {
				if(verbose) System.out.println("Time-dependent Mc = " + timeDepMcParam.getValue());
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == rmaxParam){
				
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == computeAftershockParamsButton) {
				System.out.println("Estimating aftershock model parameters...");

				tabRefreshStep = new CalcStep("Plotting", "Updating plots...", new Runnable() {

					@Override
					public void run() {
						refreshTabs(cml_num_tab_index);
					}
				}, true);

				tipStep = new CalcStep("Plotting", "Updating plots...", new Runnable() {

					@Override
					public void run() {
						printTip(4);
					}
				}, true);

				
//				CalcRunnable run = new CalcRunnable(progress, computeAftershockParamStep, computeBayesStep, postParamPlotStep, tabRefreshStep, tipStep);
				CalcRunnable run = new CalcRunnable(progress, computeBayesStep, postParamPlotStep, tabRefreshStep, tipStep);
				new Thread(run).start();
			

			} else if (param == computeAftershockForecastButton) {
				System.out.println("Computing aftershock forecast...");

				CalcStep plotStep = new CalcStep("Plotting", "Building forecast plots...",
						new Runnable() {
					@Override
					public void run() {
						plotExpectedAfershockMFDs();
						//update cumNumPlot
						plotCumulativeNum();
						setEnableParamsPostForecast(true);
					}
				}, true);

				tabSelectStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
						tabbedPane.setSelectedIndex(aftershock_expected_index);
					}

				}, true);


				tipStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
							printTip(5);
					}

				}, true);

//				CalcRunnable run = new CalcRunnable(progress, genericForecastStep, seqSpecForecastStep, bayesianForecastStep,
//						plotStep, plotTableStep, tabSelectStep, tipStep);
				CalcRunnable run = new CalcRunnable(progress, genericForecastStep, bayesianForecastStep,
						plotStep, plotTableStep, tabSelectStep, tipStep);
				new Thread(run).start();

//			} else if (param == fitMainshockSourceParam) {
//				if (fitMainshockSourceParam.getValue()) {
//					fitShakeMapSourceParam.setValue(false);
//					fitShakeMapSourceParam.getEditor().refreshParamEditor();
//				}
//
			} else if (param == fitSourceTypeParam) { 
				if (fitSourceTypeParam.getValue() == FitSourceType.CUSTOM) {
					loadSourceFromFile();
				}
				
				pgvCurves = null;
				pgaCurves = null;
				psaCurves = null;
				
				printTip(8);
				
			} else if (param == vs30Param || param == mapGridSpacingParam || param == mapScaleParam) {
				// these parameter changes are going to change the plot. Erase the old ones
				pgvCurves = null;
				pgaCurves = null;
				psaCurves = null;
				
				printTip(8);
				
			} else if (param == generateMapButton) {
				System.out.println("Computing spatial forecast map...");

				CalcStep checkValuesStep = new CalcStep("Validating", "Checking parameter values...", new Runnable() {
					
					@Override
					public void run() {
						if (mapTypeParam.getValue() == MapType.LEVEL) { 
							if (mapPOEParam.getValue() == null) 
								mapPOEParam.setValue(10);
						} else if (mapTypeParam.getValue() == MapType.PROBABILITIES){
							if (mapLevelParam.getValue() == null) 
								mapLevelParam.setValue(6);
						}
					}

				}, true);
				
				CalcStep fetchCitiesStep = new CalcStep("Building map","Loading geographical information...", new Runnable() {

					@Override
					public void run() {
//						loadCities(); //this is contained in the make map step
					}

				}, true);
				
				tipStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
						printTip(6);
					}

				}, true);

				CalcRunnable run = new CalcRunnable(progress, checkValuesStep, fetchSourceStep, fetchCitiesStep, plotMapStep, tipStep);
				new Thread(run).start();
			} else if (param == plotSpecificOnlyParam) {
				printTip(7);
				
			} else if (param == publishAdvisoryButton) {
				System.out.println("Generating advisory document...");

				CalcStep pubStep = new CalcStep("Building Document", "Building aftershock forecast summary...",
						new Runnable() {

					@Override
					public void run() {
						publishGraphicalForecast();
					}

				}, true);
			
				CalcStep quickTabRefreshStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
						refreshTabs(forecast_map_tab_index);
					}

				}, true);

				CalcRunnable run = new CalcRunnable(progress);
				if (!plotAllDurationsParam.getValue() && forecastDurationParam.getValue() >= 366d){
					CalcStep updateForecastWindowStep = new CalcStep( "Updating", "Updating foreacst duration", new Runnable() {
						@Override
						public void run() {
							plotAllDurationsParam.setValue(true);
							plotAllDurationsParam.getEditor().refreshParamEditor();
						}
					}, true);
					
					run = new CalcRunnable(progress,
							updateForecastWindowStep,
							postForecastPlotStep, 
							fetchSourceStep,
							plotMapStep,
							pubStep, 
							autoSaveStep);	
					
				} else if (forecastDurationParam.getValue() < 366d) {
					
					CalcStep updateForecastWindowStep = new CalcStep( "Updating", "Updating foreacst duration", new Runnable() {
						@Override
						public void run() {

							plotAllDurationsParam.setValue(true);
							plotAllDurationsParam.getEditor().refreshParamEditor();
							forecastDurationParam.setValue(366d);
							forecastDurationParam.getEditor().refreshParamEditor();
							updateForecastTimes();
						}
					}, true);
					
					run = new CalcRunnable(progress, 
							updateForecastWindowStep,
//							computeAftershockParamStep,
							computeBayesStep,
							postParamPlotStep,
							genericForecastStep,
//							seqSpecForecastStep,
							bayesianForecastStep,
							postForecastPlotStep,
							plotTableStep,
							fetchSourceStep,
							plotMapStep,
							quickTabRefreshStep,
							pubStep,
							autoSaveStep);
				} else
					run = new CalcRunnable(progress, fetchSourceStep, pubStep, autoSaveStep);
				
				new Thread(run).start();
			} else if (param == intensityTypeParam || param == mapTypeParam) {
//				if(verbose) System.out.println("updating map panel");
				updateMapPanel();
				
			} else if (param == mapLevelParam){
//				if(verbose) System.out.println("validating intensity value");
				
				if (intensityTypeParam.getValue() == IntensityType.MMI) {
					if (mapLevelParam.getValue() < 4)
						mapLevelParam.setValue(4d);
					else if (mapLevelParam.getValue() > 10)
						mapLevelParam.setValue(10);
				}
			} else {
				if(verbose) System.out.println("No action associated with this event.");
			}
			} catch (Exception e) {
				e.printStackTrace();
			} finally {
				setChangeListenerEnabled(true);
			}
		}
	}

	
	private void plotExpectedAfershockMFDs() {
			if (forecastMFDPane == null)
				forecastMFDPane = new JTabbedPane();
			else
				while (forecastMFDPane.getTabCount() > 0)
					forecastMFDPane.removeTabAt(0);
			
			Double minDays = forecastStartTimeParam.getValue();
			validateParameter(minDays, "start time");
			Double forecastDuration = forecastDurationParam.getValue();
			Double maxDays = minDays + forecastDuration;
			validateParameter(maxDays, "end time");
			
			double minMag;
	
			minMag = Math.min(3d, mcParam.getValue());
			
			double maxMag = 9.5d;
			
			double deltaMag = 0.1;
			int numMag = (int)((maxMag - minMag)/deltaMag + 1.5);
			
//			List<PlotElement> funcs = Lists.newArrayList();
//			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
//			
//			List<PlotElement> funcsProb = Lists.newArrayList();
//			List<PlotCurveCharacterstics> charsProb = Lists.newArrayList();
//			
			List<ETAS_AftershockModel> models = Lists.newArrayList();
			List<String> names = Lists.newArrayList();
			List<Color> colors = Lists.newArrayList();
			
//			if (seqSpecModel != null && devMode) {
//				models.add(seqSpecModel);
//				names.add("Sequence only");
//				colors.add(sequence_specific_color);
//			}
			
			if (genericModel != null && !plotSpecificOnlyParam.getValue()) {
				models.add(genericModel);
				names.add("Typical sequence");
				colors.add(generic_color);
			}
			if (bayesianModel != null) {
				models.add(bayesianModel);
				names.add("This sequence");
				colors.add(bayesian_color);
			}
					
			double[] fractiles = { 0.025, 0.500, 0.975 };

			final MinMaxAveTracker yTrack = new MinMaxAveTracker();
			final MinMaxAveTracker yTrackProb = new MinMaxAveTracker();

			int nTabs = 0;
			for (int i = 0; i < ForecastDuration.values().length; i++){
				// this logic tree decides whether to plot day.week.month.year plots or just a single plot with full duration 
				ForecastDuration  foreDur;
				String durString;
				if (plotAllDurationsParam.getValue()){
					foreDur = ForecastDuration.values()[i];
					maxDays = minDays + foreDur.duration;
					validateParameter(maxDays, "end time");
					durString = foreDur.toString();
				} else if ( i == 0 ) {
					maxDays = minDays + forecastDurationParam.getValue();
					durString = String.format("%d days", Math.round(forecastDurationParam.getValue()));
				} else {
					break;
				}
						
				
				if (maxDays <= forecastEndTimeParam.getValue()){
					
					
					List<PlotElement> funcs = Lists.newArrayList();
					List<PlotCurveCharacterstics> chars = Lists.newArrayList();
					
					List<PlotElement> funcsProb = Lists.newArrayList();
					List<PlotCurveCharacterstics> charsProb = Lists.newArrayList();
					
					
					for (int j=0; j<models.size(); j++) {
						 
						String name = names.get(j);

						if(verbose) System.out.println("Plotting forecast mfd for " + name + " - " + durString + "...");
						else System.out.format(".");

						ETAS_AftershockModel model = models.get(j);
						Color c = colors.get(j);

						EvenlyDiscretizedFunc[] fractilesFuncs = model.getCumNumMFD_FractileWithAleatoryVariability(
								fractiles, minMag, maxMag, numMag, minDays, maxDays);
						
						UncertainArbDiscDataset uncertainFunc = new UncertainArbDiscDataset(fractilesFuncs[1], fractilesFuncs[0], fractilesFuncs[2]);

						EvenlyDiscretizedFunc magPDFfunc = model.getMagnitudePDFwithAleatoryVariability(minMag, maxMag, numMag, minDays, maxDays);
						magPDFfunc.scale(100d);
						magPDFfunc.setName(name);
						funcsProb.add(magPDFfunc);
						charsProb.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));


						for (int k = 0; k < fractiles.length; k++) {
							//				double f = fractiles[k];

							EvenlyDiscretizedFunc fractile = fractilesFuncs[k];
							//				fractile.setName(name + " p"+(float)(f*100d)+"%");

							if(model!=genericModel || devMode || k != 1){

								funcs.add(fractile);
								if (k==1) {
									chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
									fractile.setName(name);
								} else {
									chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, c));
									fractile.setName("");
								}
							}
						}

						if(bayesianModel == null || devMode){
							fractilesFuncs[1].setName(name);
							funcs.add(fractilesFuncs[1]);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
						}

						uncertainFunc.setName("Range for " + name);
						funcs.add(uncertainFunc);

						PlotLineType plt = PlotLineType.SHADED_UNCERTAIN;
						Color c_trans = new Color(c.getRed(), c.getGreen(),c.getBlue(), (int) (0.3*255) );
						PlotCurveCharacterstics uncertainChars = new PlotCurveCharacterstics(plt, 1f, c_trans);
						chars.add(uncertainChars);
					}
					if(!verbose)System.out.format("\n");

					// mainshock mag and Bath's law, use evenly discr functions so that it shows up well at all zoom levels
					double mainshockMag = mainshock.getMag();
					DefaultXY_DataSet mainshockFunc = new DefaultXY_DataSet();
					DefaultXY_DataSet mainshockProbFunc = new DefaultXY_DataSet();
					mainshockProbFunc.setName("Mainshock M="+(float)mainshockMag);

					for (PlotElement elem : funcs) {
						if (elem instanceof XY_DataSet) {
							XY_DataSet xy = (XY_DataSet)elem;
							for (Point2D pt : xy)
								if (pt.getY() > 0)
									yTrack.addValue(pt.getY());
						}
					}

					for (PlotElement elem : funcsProb) {
						if (elem instanceof XY_DataSet) {
							XY_DataSet xy = (XY_DataSet)elem;
							for (Point2D pt : xy)
								if (pt.getY() > 0)
									yTrackProb.addValue(pt.getY());
						}
					}

					if (Double.isFinite(yTrack.getMin()) && Double.isFinite(yTrack.getMax()) ){
						if(D) System.out.println(yTrack);
					} else {
						yTrack.addValue(0d);
						yTrack.addValue(1d);
						if(D) System.out.println(yTrack);
					}

					yTrackProb.addValue(100d);
					if (Double.isFinite(yTrackProb.getMin()) && Double.isFinite(yTrackProb.getMax()) ){
						if(D) System.out.println(yTrackProb);
					} else {
						yTrackProb.addValue(0d);
						if(D) System.out.println(yTrackProb);
					}


					int npts = 20; 
					EvenlyDiscretizedFunc yVals = new EvenlyDiscretizedFunc(Math.min(yTrack.getMin(),1), Math.max(1e3,yTrack.getMax()), npts);
					for (int j=0; j<yVals.size(); j++) {
						double y = yVals.getX(j);
						mainshockFunc.set(mainshockMag, y);
					}

					EvenlyDiscretizedFunc yValsProb = new EvenlyDiscretizedFunc(yTrackProb.getMin(), yTrackProb.getMax(), npts);
					for (int j=0; j<yValsProb.size(); j++) {
						double y = yValsProb.getX(j);
						mainshockProbFunc.set(mainshockMag, y);
					}

					funcs.add(mainshockFunc);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
					funcsProb.add(mainshockProbFunc);
					charsProb.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));

					final PlotSpec specNum = new PlotSpec(funcs, chars, "Likely range of aftershock numbers in the next " + durString, "Magnitude", "Number of aftershocks exceeding magnitude M");
					specNum.setLegendVisible(true);

					final PlotSpec specProb = new PlotSpec(funcsProb, charsProb, "Chance of an aftershock larger than magnitude M in the next " + durString, "Magnitude", "Probability (%)");
					specProb.setLegendVisible(true);


					if (aftershockExpectedNumGraph == null){
						aftershockExpectedNumGraph = new ArrayList<GraphWidget>();
					}
					
					if(nTabs >= aftershockExpectedNumGraph.size() || aftershockExpectedNumGraph.get(nTabs) == null){
						aftershockExpectedNumGraph.add(new GraphWidget(specNum));
						GraphPanel graphPanel = aftershockExpectedNumGraph.get(nTabs).getGraphPanel();
						graphPanel.getComponent(2).setVisible(false);
					}else{
						aftershockExpectedNumGraph.get(nTabs).setPlotSpec(specNum);
					}
					aftershockExpectedNumGraph.get(nTabs).setName("Number (" + durString + ")");
					
					aftershockExpectedNumGraph.get(nTabs).setY_Log(true);
					aftershockExpectedNumGraph.get(nTabs).setY_AxisRange(new Range(1, Math.max(yTrack.getMax(),1e3)));
					aftershockExpectedNumGraph.get(nTabs).setX_AxisRange(new Range(minMag, mainshockMag+1d));
					setupGP(aftershockExpectedNumGraph.get(nTabs));

					if (aftershockProbabilityGraph == null){
						aftershockProbabilityGraph = new ArrayList<GraphWidget>();
					}	

					if(nTabs >= aftershockProbabilityGraph.size() || aftershockProbabilityGraph.get(nTabs) == null){
						aftershockProbabilityGraph.add(new GraphWidget(specProb));
						GraphPanel graphPanel = aftershockProbabilityGraph.get(nTabs).getGraphPanel();
						graphPanel.getComponent(2).setVisible(false);
					}
					else
						aftershockProbabilityGraph.get(nTabs).setPlotSpec(specProb);

					aftershockProbabilityGraph.get(nTabs).setName("Probability (" + durString + ")");
					
					aftershockProbabilityGraph.get(nTabs).setY_Log(false);
					aftershockProbabilityGraph.get(nTabs).setY_AxisRange(new Range(1e-2,100));
					aftershockProbabilityGraph.get(nTabs).setX_AxisRange(new Range(minMag, largestShock.getMag()+1d));
					setupGP(aftershockProbabilityGraph.get(nTabs));

					forecastMFDPane.insertTab("Number (" + durString + ")", null, aftershockExpectedNumGraph.get(nTabs), "Aftershock Expected Number Plot", 0);
					forecastMFDPane.insertTab("Probability (" + durString + ")", null, aftershockProbabilityGraph.get(nTabs), "Aftershock Probability Plot", 1);

					nTabs++;//add a tab
				}
			}
			if (tabbedPane.getTabCount() >= aftershock_expected_index && !tabbedPane.isEnabledAt(aftershock_expected_index)){
				tabbedPane.removeTabAt(aftershock_expected_index);
				tabbedPane.insertTab("Forecast", null, forecastMFDPane,
						"Forecast Plots", aftershock_expected_index);
			}
			else
				Preconditions.checkState(tabbedPane.getTabCount() > aftershock_expected_index, "Plots added out of order");

	}

	private void plotForecastTable(CalcProgressBar progress) {
			if(verbose) System.out.println("Building forecast tables...");
		
			if (forecastTablePane == null)
				forecastTablePane = new JTabbedPane();
			else
				while (forecastTablePane.getTabCount() > 0)
					forecastTablePane.removeTabAt(0);
			
			List<ETAS_AftershockModel> models = Lists.newArrayList();
			List<String> names = Lists.newArrayList();
	
			if (bayesianModel != null) {
				models.add(bayesianModel);
				names.add("This sequence");
			}
	
//			if (seqSpecModel != null && devMode){
//				models.add(seqSpecModel);
//				names.add("Seq. specific");
//			}
			
			if (genericModel != null) {
				models.add(genericModel);
				names.add("Typical sequence");
			}	
		
			
			GregorianCalendar eventDate = mainshock.getOriginTimeCal();
			GregorianCalendar startDate = new GregorianCalendar();
			Double minDays = forecastStartTimeParam.getValue();
			validateParameter(minDays, "start time");
			double startTime = eventDate.getTime().getTime() + minDays*ETAS_StatsCalc.MILLISEC_PER_DAY;
			startDate.setTimeInMillis((long)startTime);
			
			for (int i=0; i<models.size(); i++) {
				Stopwatch watch = Stopwatch.createStarted();
				ETAS_AftershockModel model = models.get(i);
				String name = names.get(i);
				
	//			if (progress != null)
	//				progress.updateProgress(i, models.size(), "Calculating "+name+"...");
				
				int minMag = 3;
				int maxMag = Math.max(minMag + 3, Math.min(9, (int) Math.ceil(largestShock.getMag() + 0.5)));
				int nMags = maxMag - minMag + 1;
				if(D) System.out.println("maxMag: " + maxMag + " largestShockMag: " + largestShock.getMag());
				double[] minMags = ETAS_StatsCalc.linspace(minMag, maxMag, nMags);
				
				ETAS_USGS_AftershockForecast forecast = new ETAS_USGS_AftershockForecast(model, minMags, eventDate, startDate, forecastEndTimeParam.getValue());
				JTable jTable = new JTable(forecast.getTableModel());
				jTable.getTableHeader().setFont(jTable.getTableHeader().getFont().deriveFont(Font.BOLD));
				
				forecastTablePane.addTab(name, jTable);
				
				if(D) System.out.println("Took "+watch.elapsed(TimeUnit.SECONDS)+"s to compute aftershock table for "+name);
				watch.stop();
			}
			
			if (tabbedPane.getTabCount() >= forecast_table_tab_index && !tabbedPane.isEnabledAt(forecast_table_tab_index)){
				
				
				tabbedPane.removeTabAt(forecast_table_tab_index);
				tabbedPane.insertTab("Forecast Table", null, forecastTablePane,
						"Forecast Table", forecast_table_tab_index);
			}
			else
				Preconditions.checkState(tabbedPane.getTabCount() > forecast_table_tab_index, "Plots added out of order");
		}

	private void plotRateModel2D(CalcProgressBar progress){
			
			double spacing = mapGridSpacingParam.getValue()/111.111; 	// grid spacing in degrees
			double stressDrop = 2.0; 	//MPa
			double mainshockFitDuration = dataEndTimeParam.getValue(); //days
			Double minDays = forecastStartTimeParam.getValue();
			double maxZ = 0;

			// do the full duration models
			String fitType = new String();
//			if (fitMainshockSourceParam.getValue())
			if (fitSourceTypeParam.getValue() == FitSourceType.AFTERSHOCKS)
				fitType = "aftershocks";
			else if (fitSourceTypeParam.getValue() == FitSourceType.SHAKEMAP)
				fitType = "shakemap";
			else if (fitSourceTypeParam.getValue() == FitSourceType.CUSTOM)
				fitType = "custom";
			else
				fitType = "none";
			
			// generate rate model
			GriddedGeoDataSet forecastRateModel;
			if (bayesianModel != null){
				rateModel2D = new ETAS_RateModel2D(bayesianModel);
			} else {
				rateModel2D = new ETAS_RateModel2D(genericModel);
			}

			// calculate rate for one day -- the minimum, so that we can scale up probabilities from small to large
			forecastRateModel = rateModel2D.calculateRateModel(ForecastDuration.DAY.duration, spacing, stressDrop, mainshockFitDuration, fitType, faultTrace);
			//compute rate sum for checking
			double referenceRate = 0;
			for (int j = 0; j < forecastRateModel.size(); j++){
				referenceRate += forecastRateModel.get(j);
			}
			if(D) System.out.println("rateModel sum: " + referenceRate);
			
			double value;
			boolean isMapOfProbabilities;
			if (mapTypeParam.getValue() == MapType.LEVEL) {
				value = mapPOEParam.getValue()/100; // send the % prob of exceedence as a fraction between 0 and 1
				isMapOfProbabilities = false;		// plotting levels with a given probability
			} else {
				if (intensityTypeParam.getValue() == IntensityType.PGA ||
						intensityTypeParam.getValue() == IntensityType.PSA)
					value = mapLevelParam.getValue()/100;	// send the %g as a fraction of g between 0 and 1
				else
					value = mapLevelParam.getValue();		// send MMI and PGV as is
				isMapOfProbabilities = true;		// plotting probabilities at a given level			
			}

			// calculate gmpe probabilities for one day
			GriddedGeoDataSet gmpeProbModel = getIntensityModel(forecastRateModel, value, isMapOfProbabilities, true);
			
			/*
			 * If something goews from with the intensity calculation, or if the calculation is cancelled, this is the exit point
			 */
			if (gmpeProbModel == null) return; //EXIT THE MAP COMPUTATION		
			
			
			if(D) {
			//compute rate sum for checking
				double probSum = 0;
				if (gmpeProbModel != null) {
					for (int j = 0; j < gmpeProbModel.size(); j++){
						if (!Double.isNaN(gmpeProbModel.get(j)))
							probSum += gmpeProbModel.get(j);
					}
					System.out.println("gmpeProbModel sum: " + probSum + " size: " + gmpeProbModel.size() + " max: " + gmpeProbModel.getMaxZ());
				}
			}
			
			// load up country borders
			ArrayList<XY_DataSet> countryPlot = loadCountries(gmpeProbModel.getRegion());
			
			// map city locations 
			GeoFeatureList cities = loadCities(forecastRateModel.getRegion(), 20);
			if(D) System.out.println( cities.size() + " cities returned.");
			DefaultXY_DataSet cityPlot = new DefaultXY_DataSet();
			cityPlot.setName("");
			// city names
			Font font = new Font(Font.SANS_SERIF, Font.PLAIN, 14);
			List<XYTextAnnotation> cityLabels = new ArrayList<>();
			for (GeoFeature city : cities){
				cityPlot.set(city.loc.getLongitude(), city.loc.getLatitude());
				
				XYTextAnnotation ann = new XYTextAnnotation(city.name, city.loc.getLongitude(), city.loc.getLatitude());
				ann.setTextAnchor(TextAnchor.BOTTOM_LEFT);
		        ann.setFont(font);
		        cityLabels.add(ann);
			}
			if(verbose) System.out.println( cities.size() + " cities added to plot list.");
			
			// add old hypocenter
			// scoop up the aftershock locations
			List<PlotElement> oldFuncs = new ArrayList<PlotElement>(); 
			List<PlotCurveCharacterstics> oldChars = new ArrayList<PlotCurveCharacterstics>();
			
			PlotSpec oldSpec = epicenterGraph.getPlotSpec();
			oldFuncs.addAll((List<PlotElement>) oldSpec.getPlotElems());
			oldChars.addAll(oldSpec.getChars());
			
			boolean rateModelPlotted = false;
			
			double maxContourLevel = 0;
			double minContourLevel = 0;
			double minColorLevel = 0;
			if (mapTypeParam.getValue()==MapType.PROBABILITIES) {
				minColorLevel = 0;
				minContourLevel = 0;
				maxContourLevel = 1;
			} else if (mapTypeParam.getValue()==MapType.LEVEL) {
				switch (intensityTypeParam.getValue()) {
					case MMI:
						minColorLevel = 4;
						minContourLevel = 4;
						maxContourLevel = 5;
						break;
					case PGA:
						minColorLevel = 0;
						minContourLevel = 3;
						maxContourLevel = 4;
						break;
					case PGV:
						minColorLevel = 0;
						minContourLevel = 2;
						maxContourLevel = 3;
						break;
					case PSA:
						minColorLevel = 0;
						minContourLevel = 3;
						maxContourLevel = 4;
						break;
					default:
						minColorLevel = 0;
						minContourLevel = 0;
						maxContourLevel = 1;
						break;
					}
			}
				
			// moved here so that canceling map generation will still leave you with the previous maps.
			if (forecastMapPane == null)
				forecastMapPane = new JTabbedPane();
			else
				while (forecastMapPane.getTabCount() > 0)
					forecastMapPane.removeTabAt(0);
	
			// clear out (or initialize) the contourModel
			contourList = new ArrayList<ContourModel>();
			

			// for each duration (day/week/month/year or just what's specified) 
			for (int i = 0; i < ForecastDuration.values().length; i++){
				// clone the gridded rate model
//				GriddedGeoDataSet newForecastRateModel = forecastRateModel.copy();
				GriddedGeoDataSet newForecastRateModel = cloneGriddedGeoDataSet(forecastRateModel);
				
				// this logic tree decides whether to plot day.week.month.year plots or just a single plot with full duration 
				ForecastDuration foreDur;
				String durString;
				double plotDur, maxDays;
				if (plotAllDurationsParam.getValue()){
					foreDur = ForecastDuration.values()[i];
					plotDur = foreDur.duration;
					maxDays = minDays + plotDur;
					//					validateParameter(maxDays, "end time");
					durString = foreDur.toString();
				} else if ( i == 0 ) {
					plotDur = forecastDurationParam.getValue();
					maxDays = minDays + plotDur;
					durString = String.format("%d days", Math.round(forecastDurationParam.getValue()));
				} else {
					break;
				}
				

				if (maxDays <= forecastEndTimeParam.getValue()){
					if(verbose) System.out.println("Computing map for " + durString);
					
					//rescale the rate model (for speed) to get a new probability model
					double targetRate;
					if (bayesianModel == null){
						targetRate = genericModel.getCumNumFractileWithAleatory(new double[]{0.5}, genericModel.magComplete, minDays, maxDays)[0];
					} else {
						targetRate =  bayesianModel.getCumNumFractileWithAleatory(new double[]{0.5}, bayesianModel.magComplete, minDays, maxDays)[0];
					}
					if(D) System.out.println("rateModel & targetRate (" + i + "): " + referenceRate + " " + targetRate);
					newForecastRateModel.scale(targetRate/referenceRate);
					if(D) System.out.println("newForecastRateModel sum: " + newForecastRateModel.getSumZ());
					
					GriddedGeoDataSet new_gmpeProbModel;

					/*
					 * mess things up:
					 */
//					pgvCurves = null;
//					pgaCurves = null;
//					psaCurves = null;
					
					if (mapTypeParam.getValue() == MapType.LEVEL) //if we're plotting levels with a given prob, start from scratch
//						new_gmpeProbModel = getIntensityModel(newForecastRateModel, value, false);
						new_gmpeProbModel = scaleLevelModel(gmpeProbModel.getRegion(), referenceRate, targetRate, value);
					else { //if we're plotting probs at a given level, just scale it up
//						new_gmpeProbModel = getIntensityModel(newForecastRateModel, value, true);
						new_gmpeProbModel = scaleProbabilityModel(gmpeProbModel, referenceRate, targetRate);
					}
					
					// plot the rate only the first time around (should be the longest duration, and all the others are identical except for scaling)
					if (!rateModelPlotted) {
						rateModelPlotted = true;
						
//						GriddedGeoDataSet scaledRateModel = forecastRateModel.copy();
						GriddedGeoDataSet scaledRateModel = cloneGriddedGeoDataSet(forecastRateModel);
						if (D) System.out.println("forecastRateModel sum: " + forecastRateModel.getSumZ());
						
						
						double magRef = 3.5;
						double magScaleFactor = Math.pow(10, -bParam.getValue()*(magRef - mcParam.getValue()));
						double gridRef = 1;
						double gridScaleFactor = Math.pow(gridRef/mapGridSpacingParam.getValue(),2);
						scaledRateModel.scale(magScaleFactor*gridScaleFactor);
						
						if (D) System.out.println("forecastRateModel sum: " + forecastRateModel.getSumZ() + " after scaling clone");
						if (D) System.out.println("scaledRateModel sum: " + scaledRateModel.getSumZ());
						
						maxZ = scaledRateModel.getMaxZ();
						CPT cpt = getProbCPT().rescale(0, maxZ);
						
						List<PolyLine> rateContours = ETAS_RateModel2D.getContours(scaledRateModel, 20);
						List<PlotElement> funcs = Lists.newArrayList();
						List<PlotCurveCharacterstics> chars = Lists.newArrayList();
						
						// add epicenters to map
						for (int j = 0; j < oldFuncs.size(); j++ ){
							if (oldChars.get(j).getSymbol() == PlotSymbol.CIRCLE){
								funcs.add(oldFuncs.get(j));
								chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, (float) (oldChars.get(j).getSymbolWidth()/10), Color.GRAY));
							}
						}
						
						// add contours to map
						for (int n = 0; n < rateContours.size(); n++){
							XY_DataSet contour = new DefaultXY_DataSet();
							for (PointD pt : rateContours.get(n).PointList)
								contour.set(pt.X,pt.Y);

							funcs.add(contour);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, cpt.getColor((float) rateContours.get(n).Value)));
						}
//						Collections.reverse(funcs);
//						Collections.reverse(chars);

						// add cities to map
						funcs.add(cityPlot);
						chars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 2, Color.BLACK));
						
						// add countries to map
						for (XY_DataSet elem : countryPlot) {
							funcs.add(elem);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1, Color.GRAY));
						}
						
						PlotSpec spec = new PlotSpec(funcs, chars, "Past Aftershocks (dots) and expected number of M≥" + String.format("%2.1f", magRef)  + " aftershocks\n per square km in the next day (contours)", "Longitude", "Latitude");
//						PlotSpec spec = new PlotSpec(funcs, chars, "Probable location of the next aftershock", "Longitude", "Latitude");
						
						spec.setPlotAnnotations(cityLabels);
						
						PaintScaleLegend subtitle = XYZGraphPanel.getLegendForCPT(cpt, "Number", axisLabelFontSize, tickLabelFontSize,
					             0, RectangleEdge.RIGHT);
						if (subtitle != null)
						spec.addSubtitle(subtitle);

						// try (and fail) to disable auto-include zero 
						GraphWidget forecastMapGraph = new GraphWidget(spec);
						((NumberAxis) forecastMapGraph.getGraphPanel().getXAxis()).setAutoRangeIncludesZero(false);;
						((NumberAxis) forecastMapGraph.getGraphPanel().getYAxis()).setAutoRangeIncludesZero(false);;
						
						// clean up the graph panel
						Component buttonPanel = forecastMapGraph.getButtonControlPanel();
						buttonPanel.setVisible(false);
						GraphPanel graphPanel = forecastMapGraph.getGraphPanel();
						graphPanel.getComponent(2).setVisible(false);
						setupGP(forecastMapGraph);

						// set the axis ranges
						double regBuff = 0.05;
						Region region = newForecastRateModel.getRegion();
						forecastMapGraph.setAxisRange(region.getMinLon()-regBuff, region.getMaxLon()+regBuff,
								region.getMinLat()-regBuff, region.getMaxLat()+regBuff);
					

						forecastMapGraph.setName("Rate Map");
						forecastMapPane.addTab("Rate Map", null, forecastMapGraph);
					}
					
					if (new_gmpeProbModel != null && new_gmpeProbModel.size() > 0) {
						// draw probGraph. If the map is probabilities, or percent g, scale it up by a factor of 100 for plotting %
						if (mapTypeParam.getValue() == MapType.PROBABILITIES || intensityTypeParam.getValue() == IntensityType.PGA || 
								intensityTypeParam.getValue() == IntensityType.PSA)
							new_gmpeProbModel.scale(100d);

						// update the max contour level for setting the colorscale
						if (new_gmpeProbModel.getMaxZ() > maxContourLevel) 
							maxContourLevel = new_gmpeProbModel.getMaxZ();

						// set up the contours, with colors
						CPT cpt2;
						double[] contourLevels;
						if (mapTypeParam.getValue()==MapType.PROBABILITIES) {
							double cptMax;
							
							if (maxContourLevel > 20)
								cptMax = Math.min(100, (Math.ceil(maxContourLevel/20))*20);
							else if (maxContourLevel > 2)
								cptMax = Math.min(100, (Math.ceil(maxContourLevel/2))*2);
							else
								cptMax = 2;
															
							cpt2 = getProbCPT().rescale(minColorLevel, cptMax);
							contourLevels = ETAS_StatsCalc.linspace(minContourLevel,Math.max(new_gmpeProbModel.getMaxZ(),cptMax),21);	//specify contourLevels directly
						} else {
							cpt2 = getProbCPT().rescale(minColorLevel, Math.max(new_gmpeProbModel.getMaxZ(),maxContourLevel));
							contourLevels = ETAS_StatsCalc.linspace(minContourLevel,Math.max(new_gmpeProbModel.getMaxZ(),maxContourLevel),21);	//specify contourLevels directly
						}
						List<PolyLine> contours = ETAS_RateModel2D.getContours(new_gmpeProbModel, contourLevels);
						
						
						// spit out the contour levels used
						double maxPlottedContour = 0;
							StringBuilder outString = new StringBuilder();
							for (PolyLine line : contours) {
								outString.append(String.format("%2.2f ", line.Value));
								if (line.Value > maxPlottedContour) maxPlottedContour = line.Value;
							}
							if(verbose){ System.out.println("contour levels: " + outString);
						}

						// set up the map title
						String title = "";
						String units = "";
						String contourUnits = "";
						if (mapTypeParam.getValue() == MapType.LEVEL) {
							String level;
							if (intensityTypeParam.getValue() == IntensityType.MMI) {
								if (maxPlottedContour > 0)
									level = decToRoman((int) (Math.floor(maxPlottedContour)));
								else
									level = "<" + decToRoman((int) minContourLevel);
							} else {
								if (maxPlottedContour > 0)
									level = String.format("%.0f", maxPlottedContour);
								else
									level = String.format("<%.0f", minContourLevel);
							}
							
							
							if (intensityTypeParam.getValue() == IntensityType.PGA ||
									intensityTypeParam.getValue() == IntensityType.PSA) 
								units = "%g";
							else if (intensityTypeParam.getValue() == IntensityType.PGV)
								units = "cm/s";
							else
								units = "";
						
							String intensityType = intensityTypeParam.getValue().getAbbreviation();
							String probability = String.format("%.0f", mapPOEParam.getValue());
							title = intensityType + " with " + probability + "% chance of being exceeded in the next " + durString + " (max: " + level + units + ")";
							units = "%";
						} else if (mapTypeParam.getValue() == MapType.PROBABILITIES) {
							String intensityType = intensityTypeParam.getValue().getAbbreviation();
							String level;
							if (intensityTypeParam.getValue() == IntensityType.MMI)
								level = decToRoman((int) (Math.floor(mapLevelParam.getValue() + 0.01)));
							else
								level = String.format("%.0f", mapLevelParam.getValue());
							
							if (intensityTypeParam.getValue() == IntensityType.PGA ||
									intensityTypeParam.getValue() == IntensityType.PSA) 
								units = "%g";
							else if (intensityTypeParam.getValue() == IntensityType.PGV)
								units = "cm/s";
							else
								units = "";
							title = "Probability of exceeding " + intensityType + " level " + level + units + " in the next " + durString + " (max: " + ((maxPlottedContour < 1)?"<1":String.format("%.0f", maxPlottedContour)) + " " +"%)";
							contourUnits = "%";
						}
						contourList.add(new ContourModel(contours, title, contourUnits, cpt2));

						//copy epicenter graph and add to map
						List<PlotElement> funcs = Lists.newArrayList();
						List<PlotCurveCharacterstics> chars = Lists.newArrayList();
						for (int j = 0; j < oldFuncs.size(); j++ ){
							if (oldChars.get(j).getSymbol() == PlotSymbol.CIRCLE){
								funcs.add(oldFuncs.get(j));
								chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE,  (float) (oldChars.get(j).getSymbolWidth()/10), Color.GRAY));
							}
						}

						// add contours to map 
						for (int n = 0; n < contours.size(); n++){
							XY_DataSet contour = new DefaultXY_DataSet();
							for (PointD pt : contours.get(n).PointList)
								contour.set(pt.X,pt.Y);

							funcs.add(contour);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, cpt2.getColor((float) contours.get(n).Value)));
						}
						//					Collections.reverse(funcs);
						//					Collections.reverse(chars);

						// add the cities to map
						// re-map city locations for larger region. 
						cities = loadCities(new_gmpeProbModel.getRegion(), 20);
						if(D) System.out.println( cities.size() + " cities returned.");
						cityPlot = new DefaultXY_DataSet();
						cityPlot.setName("");
						// city names
						font = new Font(Font.SANS_SERIF, Font.PLAIN, 14);
						cityLabels = new ArrayList<>();
						for (GeoFeature city : cities){
							cityPlot.set(city.loc.getLongitude(), city.loc.getLatitude());
							
							XYTextAnnotation ann = new XYTextAnnotation(city.name, city.loc.getLongitude(), city.loc.getLatitude());
							ann.setTextAnchor(TextAnchor.BOTTOM_LEFT);
					        ann.setFont(font);
					        cityLabels.add(ann);
						}
						if(D) System.out.println( cities.size() + " cities added to plot list.");
						funcs.add(cityPlot);
						chars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 2, Color.BLACK));

						// add countries to map
						for (XY_DataSet elem : countryPlot) {
							funcs.add(elem);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1, Color.GRAY));
						}						
						PlotSpec spec = new PlotSpec(funcs, chars, title, "Longitude", "Latitude");
						spec.setPlotAnnotations(cityLabels);

						// set up the color axis plot label
						String cptAxisLabel = "";
						if (mapTypeParam.getValue() == MapType.LEVEL) {
							String cptUnits;
							if (intensityTypeParam.getValue() == IntensityType.PGA ||
									intensityTypeParam.getValue() == IntensityType.PSA) 
								cptUnits = "Peak Acceleration (%g)";
							else if (intensityTypeParam.getValue() == IntensityType.PGV)
								cptUnits = "Peak Velocity (cm/s)";
							else
								cptUnits = "MMI";

							cptAxisLabel = cptUnits;
						} else if (mapTypeParam.getValue() == MapType.PROBABILITIES) {
							String cptUnits = "Probability (%)";
							cptAxisLabel = cptUnits;
						}

						PaintScaleLegend subtitle = XYZGraphPanel.getLegendForCPT(cpt2, cptAxisLabel, axisLabelFontSize, tickLabelFontSize,
								0, RectangleEdge.RIGHT);
						if (subtitle != null)
							spec.addSubtitle(subtitle);

						GraphWidget forecastGraph = new GraphWidget(spec);
						
						// try to remove the auto-include-zero thing. Doesn't work though.
						forecastGraph.getGraphPanel().getXAxis().setDefaultAutoRange(forecastGraph.getGraphPanel().getXAxis().getRange());
						forecastGraph.getGraphPanel().getYAxis().setDefaultAutoRange(forecastGraph.getGraphPanel().getYAxis().getRange());
						
						// clean up the graph panel
						Component buttonPanel = forecastGraph.getButtonControlPanel();
						buttonPanel.setVisible(false); // get rid of the button panel
						GraphPanel graphPanel = forecastGraph.getGraphPanel();
						graphPanel.getComponent(2).setVisible(false);	//get rid of the other stuff too

						// set axis ranges
						double regBuff = 0.05;
						Region region = new_gmpeProbModel.getRegion();
						forecastGraph.setAxisRange(region.getMinLon()-regBuff, region.getMaxLon()+regBuff,
								region.getMinLat()-regBuff, region.getMaxLat()+regBuff);

						setupGP(forecastGraph);
						
						forecastGraph.setName("Shaking (" + durString + ")");
						forecastMapPane.addTab("Shaking (" + durString + ")", null, forecastGraph);
					}
						
					if (tabbedPane.getTabCount() >= forecast_map_tab_index && !tabbedPane.isEnabledAt(forecast_map_tab_index)){

						tabbedPane.removeTabAt(forecast_map_tab_index);
						tabbedPane.insertTab("Forecast Maps", null, forecastMapPane,
								"Forcast Maps (rate, intensity)", forecast_map_tab_index);
					}
					else
						Preconditions.checkState(tabbedPane.getTabCount() > forecast_map_tab_index, "Plots added out of order");
				}
			}
	}

	GriddedGeoDataSet getMMI_IntensityModel(GriddedGeoDataSet rateModel, double value, boolean isMapOfProbabilities, boolean prompt) {
		
		if (isMapOfProbabilities) {
			// Plotting probabilities. Value is MMI level, figure out what PGV and PGA levels to be looking for. 
		 
			double pga = Wald_MMI_Calc.getPGA(value);
			double pgv = Wald_MMI_Calc.getPGV(value);
			double weightPGV = Wald_MMI_Calc.getWeightVMMI(value);
			
			if (D) System.out.println("MMI level: " + value + " PGA level: " + pga + " PGV level: " + pgv + " PGV weight: " + weightPGV);
			// first get the PGA model
			intensityTypeParam.setValue(IntensityType.PGA);
			GriddedGeoDataSet pgaMap = getIntensityModel(rateModel, pga, isMapOfProbabilities, prompt);
			if (pgaMap == null)
				return null;
			
			// then get the PGV model
			intensityTypeParam.setValue(IntensityType.PGV);
			GriddedGeoDataSet pgvMap = getIntensityModel(rateModel, pgv, isMapOfProbabilities, false);
			if (pgvMap == null)
				return null;
			
			// put things back the way they were
//			mapLevelParam.setValue(level);
			intensityTypeParam.setValue(IntensityType.MMI);

			GriddedGeoDataSet mmiMap = new GriddedGeoDataSet(pgaMap.getRegion(), pgaMap.isLatitudeX());
			for (int index=0; index<mmiMap.size(); index++) {
				double pgaProb = pgaMap.get(index);
				double pgvProb = pgvMap.get(index);
				double mmiProb = calcMMI_POE(weightPGV, pgaProb, pgvProb);
				mmiMap.set(index, mmiProb);
			}
			return mmiMap;
		} else {
			// 	plotting levels. Value is POE just compute the PGV level and PGA level independently, and combine
			
 			// first get the PGA model
			if(D) System.out.println("Calculating PGA map with " + value + " probability of being exceeded...");
			intensityTypeParam.setValue(IntensityType.PGA);
			GriddedGeoDataSet pgaMap = getIntensityModel(rateModel, value, isMapOfProbabilities, prompt);

			// then get the PGV model
			if(D) System.out.println("Calculating PGV map with " + value + " probability of being exceeded...");
			intensityTypeParam.setValue(IntensityType.PGV);
			GriddedGeoDataSet pgvMap = getIntensityModel(rateModel, value, isMapOfProbabilities, false);
			
			// put things back the way they were
			intensityTypeParam.setValue(IntensityType.MMI);

			GriddedGeoDataSet mmiMap = new GriddedGeoDataSet(pgaMap.getRegion(), pgaMap.isLatitudeX());
			if(D) System.out.println("Combining PGA and PGV to get MMI...");
			for (int index=0; index<mmiMap.size(); index++) {
				double pga = pgaMap.get(index);
				if (!Doubles.isFinite(pga))
					pga = 0;
				double pgv = pgvMap.get(index);
				if (!Doubles.isFinite(pgv))
					pgv = 0;
				double mmi = Wald_MMI_Calc.getMMI(pga, pgv);
				if(!Doubles.isFinite(mmi))
					if(D) System.err.println(String.format("Bad MMI=%s for PGA=%s, PGV=%s. MMI map may have errors", mmi, pga, pgv));
				
				if (mmi == 1d)
					// will speed things up
					mmi = Double.NaN;

				mmiMap.set(index, mmi);
			}
			return mmiMap;
			
		}
		
	}
	
	private double calcMMI_POE(double weightPGV, double pgaProb, double pgvProb) {
		if (!Doubles.isFinite(pgaProb))
			pgaProb = 0;
		if(! ((float)pgaProb <= 1f))
			if(D) System.err.println(String.format("Bad PGA prob: %s", pgaProb));
		
		if (!Doubles.isFinite(pgvProb))
			pgvProb = 0;
		if(! ((float)pgvProb <= 1f))
			if(D) System.err.println(String.format("Bad PGV prob: %s", pgvProb));
		
		double mmiProb = weightPGV*pgvProb + (1d-weightPGV)*pgaProb;
		if(!(mmiProb >= 0d && mmiProb <= 1d) ) 
			if(D) System.err.println(String.format("Bad MMI probability (%s) for pgaProb= %s and pgvProb=%s",
				mmiProb, pgaProb, pgvProb));
		
		if (mmiProb == 0d)
			// will speed things up
			mmiProb = Double.NaN;
		
		return mmiProb;
	}
	
	
	
	GriddedGeoDataSet getIntensityModel(GriddedGeoDataSet rateModel, double value, boolean isMapOfProbabilities) {
		return getIntensityModel(rateModel, value, isMapOfProbabilities, false);
	}
	
	GriddedGeoDataSet getIntensityModel(GriddedGeoDataSet rateModel, double value, boolean isMapOfProbabilities, boolean prompt) {
		if(verbose) System.out.println("Computing " + (isMapOfProbabilities?"probability":"level") + " map for type: " +
					intensityTypeParam.getValue().getAbbreviation() + " at level: " + value);
		
		//MMI is sort of implemented, by running PGA and PGV both
		if (intensityTypeParam.getValue() == IntensityType.MMI) {
			// short circuit to the MMI method
			return getMMI_IntensityModel(rateModel, value, isMapOfProbabilities, prompt);			
		}
		
		// first order of business is to correct the rate model to the "allowable" magnitude range. 
		double mc = mcParam.getValue();
		double refMag = 4.0;
		// clone the gridded rate model
//		GriddedGeoDataSet newRateModel = rateModel.copy();
		GriddedGeoDataSet newRateModel = cloneGriddedGeoDataSet(rateModel);
		
		if(D) System.out.println("mc: " + mc + " mref: " + refMag);
		if (mc < refMag)
			newRateModel.scale(Math.pow(10, -bParam.getValue()*(refMag - mc)));
		else
			refMag = mc;

		/*
		 * debugging. choose the top one for the real version
		 */
		// set up an expanded region, based on user-specified range factor
		double currentRadius = 111.111*(rateModel.getRegion().getMaxLat() - rateModel.getRegion().getMinLat())/2d;
		double calcSpacing = rateModel.getRegion().getSpacing()*mapScaleParam.getValue()*2;		
		SphRegionCircle newRegion = new SphRegionCircle(new SphLatLon(getCentroid()), currentRadius*mapScaleParam.getValue());
		GriddedRegion calcRegion = new GriddedRegion(new Region(new Location(newRegion.getMaxLat(), newRegion.getMaxLon()),
				new Location(newRegion.getMinLat(), newRegion.getMinLon())), calcSpacing, null);
		
		if(D) System.out.println("currentRadius: " + currentRadius + " calcSpacing: " + calcSpacing);
		//this is what is in the old code
//		double calcSpacing = rateModel.getRegion().getSpacing()*2;
//		GriddedRegion calcRegion = new GriddedRegion(new Region(new Location(rateModel.getMaxLat(), rateModel.getMaxLon()),
//				new Location(rateModel.getMinLat(), rateModel.getMinLon())), calcSpacing, null);
		/*
		 * 
		 */
		
		ScalarIMR gmpe = AttenRelRef.BSSA_2014.instance(null); //consider which attenuation relationship to use...
//		ScalarIMR gmpe = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null); //consider which attenuation relationship to use...
		gmpe.setParamDefaults();

		if (intensityTypeParam.getValue().equals(IntensityType.PGA)) {
			gmpe.setIntensityMeasure(PGA_Param.NAME);
		} else if (intensityTypeParam.getValue().equals(IntensityType.PGV)) { 
			gmpe.setIntensityMeasure(PGV_Param.NAME);
		} else if (intensityTypeParam.getValue().equals(IntensityType.PSA)) {
			gmpe.setIntensityMeasure(SA_Param.NAME);
			SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), 1d);
		} else {
			// no map, return empty data set
			return new GriddedGeoDataSet(calcRegion, false);
		}
		
		
		// Vs30 provider, or null for no Vs30
		WaldAllenGlobalVs30 vs30Provider;
		if (!vs30Param.getValue()) {
			vs30Provider = null;
		} else {
			try {
				vs30Provider = new WaldAllenGlobalVs30();
				TectonicRegime regime = tectonicRegimeParam.getValue();
				if (regime == TectonicRegime.SCR_ABVSLAB ||
						regime == TectonicRegime.SCR_GENERIC ||
						regime == TectonicRegime.SOR_ABVSLAB ||
						regime == TectonicRegime.SOR_GENERIC) {
					vs30Provider.setStableCoefficients();
				} else {
					vs30Provider.setActiveCoefficients();
				}
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				System.err.println("Could not retrieve global DEM for estimating site effects. Proceeding without site effects.");
				e1.printStackTrace();
				return null;
			}
		}

		// this is where you can code the tectonic style/dominant focal mechanism -- hardcoded for now.
		Map<FocalMech, Double> mechWts = new HashMap<>();
		mechWts.put(FocalMech.STRIKE_SLIP, 0.5);
		mechWts.put(FocalMech.NORMAL, 0.25);
		mechWts.put(FocalMech.REVERSE, 0.25);

	
		
//		GriddedRegion calcRegion = new GriddedRegion(new Region(new Location(rateModel.getMaxLat(), rateModel.getMaxLon()),
//				new Location(rateModel.getMinLat(), rateModel.getMinLon())), calcSpacing, null);

		
		
		double maxSourceDist = 200d;
		
		for (Location calcLoc : calcRegion.getBorder())
			for (Location rateLoc : rateModel.getRegion().getBorder())
				maxSourceDist = Math.max(maxSourceDist, LocationUtils.horzDistanceFast(calcLoc, rateLoc));
		if(D) System.out.println("Maximum source distance: "+(float)maxSourceDist);

		DiscretizedFunc[] curves;
		GriddedGeoDataSet map;
		
		double minMag = refMag;
		double maxMag = Math.min(genericParams.get_maxMag(), 9);
		
		try {
			if (intensityTypeParam.getValue()==IntensityType.PGV) {
				if (pgvCurves == null) {
					if(D) System.out.println("Calculating pgv curves at each grid point");
					curves = ETAS_ShakingForecastCalc.calcForecast(calcRegion, newRateModel, minMag, maxMag, bParam.getValue(), gmpe, mechWts,
							maxSourceDist, vs30Provider, prompt);
					if (curves != null) pgvCurves = curves.clone();
				} else {
					if(D) System.out.println("reloading pgv curves");
					curves = pgvCurves;
				}
			} else if (intensityTypeParam.getValue()==IntensityType.PGA) {
				if (pgaCurves == null) {
					if(D) System.out.println("Calculating pga curves at each grid point");
					curves = ETAS_ShakingForecastCalc.calcForecast(calcRegion, newRateModel, minMag, maxMag, bParam.getValue(), gmpe, mechWts,
							maxSourceDist, vs30Provider, prompt);
					if (curves != null) pgaCurves = curves.clone();
				} else {
					if(D) System.out.println("reloading pga curves");
					curves = pgaCurves;
				}
			} else if (intensityTypeParam.getValue()==IntensityType.PSA) {
				if (psaCurves == null) {
					if(D) System.out.println("Calculating pga curves at each grid point");
					curves = ETAS_ShakingForecastCalc.calcForecast(calcRegion, newRateModel, minMag, maxMag, bParam.getValue(), gmpe, mechWts,
							maxSourceDist, vs30Provider, prompt);
					if (curves != null) psaCurves = curves.clone();
				} else {
					if(D) System.out.println("reloading psa curves");
					curves = psaCurves;
				}
			} else {
				curves = null;
			}
			
			if (curves == null) {
				if(verbose) System.out.println("Cancelling map generation...");
					return null;
			} else {
				if(D) System.out.println("Extracting map for value: " + value + ". Plotting probabilities? " + isMapOfProbabilities);
				map = ETAS_ShakingForecastCalc.extractMap(calcRegion, curves, isMapOfProbabilities, value); 
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		
		// filter the map for NaN values
		int nanSum = 0;
		for (int i =0; i < map.size(); i++) {
			if(Double.isNaN(map.get(i))) {
				map.set(i, 0d);
				nanSum++;
			}
		}
		
		if(D) System.out.println(nanSum + " NaN values found in probability model");
		return map;

	}

	private GriddedGeoDataSet scaleMMILevelModel(GriddedRegion region, double referenceRate, double targetRate, double value) {
		// 	plotting levels. Value is POE just compute the PGV level and PGA level independently, and combine
		
		// first get the PGA model
		if(D) System.out.println("Calculating PGA map with " + value + " probability of being exceeded");
		intensityTypeParam.setValue(IntensityType.PGA);
		GriddedGeoDataSet pgaMap = scaleLevelModel(region, referenceRate, targetRate, value);

		// then get the PGV model
		if(D) System.out.println("Calculating PGV map with " + value + " probability of being exceeded");
		intensityTypeParam.setValue(IntensityType.PGV);
		GriddedGeoDataSet pgvMap = scaleLevelModel(region, referenceRate, targetRate, value);
		
		// put things back the way they were
		intensityTypeParam.setValue(IntensityType.MMI);

		GriddedGeoDataSet mmiMap = new GriddedGeoDataSet(pgaMap.getRegion(), pgaMap.isLatitudeX());
		if(D) System.out.println("Combining PGA and PGV to get MMI");
		for (int index=0; index<mmiMap.size(); index++) {
			double pga = pgaMap.get(index);
			if (!Doubles.isFinite(pga))
				pga = 0;
			double pgv = pgvMap.get(index);
			if (!Doubles.isFinite(pgv))
				pgv = 0;
			double mmi = Wald_MMI_Calc.getMMI(pga, pgv);
			if(!Doubles.isFinite(mmi))
				if(D) System.err.println(String.format("Bad MMI=%s for PGA=%s, PGV=%s. MMI map may have errors", mmi, pga, pgv));
			
			if (mmi == 1d)
				// will speed things up
				mmi = Double.NaN;

			mmiMap.set(index, mmi);
		}
		return mmiMap;
			
	}
	
	private GriddedGeoDataSet scaleLevelModel(GriddedRegion region, double referenceRate, double targetRate, double value) {
		
		DiscretizedFunc[] curves;
		if (intensityTypeParam.getValue()==IntensityType.PGV) 
			curves = pgvCurves;
		else if (intensityTypeParam.getValue()==IntensityType.PGA)
			curves = pgaCurves;
		else if (intensityTypeParam.getValue()==IntensityType.PSA)
			curves = psaCurves;
		else if (intensityTypeParam.getValue()==IntensityType.MMI)
			return scaleMMILevelModel(region, referenceRate, targetRate, value);
		else if (intensityTypeParam.getValue()==IntensityType.NONE)
			return null;
		else
			return null;
		
		if(D) System.out.println("referenceRate: " + referenceRate + " targetRate: " + targetRate);
		DiscretizedFunc[] newCurves = new DiscretizedFunc[curves.length];
		//scale the curves up
		for (int i = 0; i < curves.length; i++) {
			//convert to rate
			DiscretizedFunc rateCurve = poissProb2Rate(curves[i].deepClone());
			rateCurve.scale(targetRate/referenceRate);
			newCurves[i] = rate2PoissProb(rateCurve);
		}
		
		GriddedGeoDataSet map = ETAS_ShakingForecastCalc.extractMap(region, newCurves, false, value);
		return map;
	}
	
	private DiscretizedFunc rate2PoissProb(DiscretizedFunc rateModel){
		double prob;
		for (int i = 0; i < rateModel.size(); i++){
			prob = 1 - Math.exp(-rateModel.getY(i));
			rateModel.set(i, Math.min(prob, 0.999));
		}
		
		return rateModel;
	}
	
	private DiscretizedFunc poissProb2Rate(DiscretizedFunc probModel){
		double rate;
		for (int i = 0; i < probModel.size(); i++){
			rate = -Math.log(1-probModel.getY(i));
			probModel.set(i, rate);
		}
		
		return probModel;
	}
	
	
	private GriddedGeoDataSet scaleProbabilityModel(GriddedGeoDataSet gmpeProbModel, double referenceRate, double targetRate) {		
		 if (intensityTypeParam.getValue()==IntensityType.NONE)
			return null;
		 
		// the gmpeRateModel is the prob model converted to Poisson rate 
//		GriddedGeoDataSet gmpeRateModel = gmpeProbModel.copy();
		GriddedGeoDataSet gmpeRateModel = cloneGriddedGeoDataSet(gmpeProbModel);

		poissProb2Rate(gmpeRateModel); //convert the probabilities to a rate, for scaling.

//		//compute rate sum for reference model
//		double gmpeRateSum = 0;
//		for (int j = 0; j < gmpeRateModel.size(); j++){
//			if (!Double.isNaN(gmpeRateModel.get(j)))
//				gmpeRateSum += gmpeRateModel.get(j);
//		}
		
//		if(D) System.out.println("gmpeRateModel sum: " + gmpeRateSum + " max: " + gmpeProbModel.getMaxZ());

//		gmpeRateModel.scale(targetRate/gmpeRateSum);
		gmpeRateModel.scale(targetRate/referenceRate);
		GriddedGeoDataSet new_gmpeProbModel = rate2PoissProb(gmpeRateModel);

		return new_gmpeProbModel;
	}

	private GriddedGeoDataSet rate2PoissProb(GriddedGeoDataSet rateModel){
		double prob;
		for (int i = 0; i < rateModel.size(); i++){
			prob = 1 - Math.exp(-rateModel.get(i));
			rateModel.set(i, Math.min(prob,0.999));
		}
		
		return rateModel;
	}
		
	private GriddedGeoDataSet poissProb2Rate(GriddedGeoDataSet probModel){
		double rate;
		for (int i = 0; i < probModel.size(); i++){
			rate = -Math.log(1-probModel.get(i));
			probModel.set(i, rate);
		}
		
		return probModel;
	}
	
	private String decToRoman(int mmi){
		String str;
		switch(mmi) {
			case 0:
				str = "I";
				break;
			case 1:
				str = "I";
				break;
			case 2:
				str = "II";
				break;
			case 3:
				str = "III";
				break;
			case 4:
				str = "IV";
				break;
			case 5:
				str = "V";
				break;
			case 6:
				str = "VI";
				break;
			case 7:
				str = "VII";
				break;
			case 8:
				str = "VIII";
				break;
			case 9:
				str = "IX";
				break;
			case 10:
				str = "X";
				break;
			default:
				str = "X";
				break;
		}
		return str;
	}
	
	private void publishGraphicalForecast(){
		// use a popup or dropdown menu to set graphical forecast fields
		
		//this is cribbed from the previous method
		GregorianCalendar eventDate = mainshock.getOriginTimeCal();
		GregorianCalendar startDate = new GregorianCalendar();
		Double minDays = forecastStartTimeParam.getValue();
		validateParameter(minDays, "start time");
		double startTime = eventDate.getTime().getTime() + minDays*ETAS_StatsCalc.MILLISEC_PER_DAY;
		startDate.setTimeInMillis((long)startTime);
				
		// write HTML document summary 
		File outFile;
		workingDirChooser = new JFileChooser(workingDir.getPath());
		workingDirChooser.setSelectedFile(new File("Advisory.html"));
		FileNameExtensionFilter htmlFilter = new FileNameExtensionFilter("html files (*.html)", "html");
		workingDirChooser.setFileFilter(htmlFilter);
		
		int ret = workingDirChooser.showSaveDialog(this);
		
		if (ret == JFileChooser.APPROVE_OPTION) {
			outFile = workingDirChooser.getSelectedFile();
			workingDir = outFile.getParentFile();
			workingDirChooser.setCurrentDirectory(workingDir);
			workingDirChooser.changeToParentDirectory();
			
			ETAS_AftershockModel model;
			if (bayesianModel != null)
				model = bayesianModel;
			else 
				model = genericModel;
			
			if(verbose) System.out.println("Publishing forecast...");
			
			GraphicalForecast graphForecast = new GraphicalForecast(outFile, model, eventDate, startDate);
			graphForecast.setShakeMapURL(shakeMapURL);
			graphForecast.constructForecast();
			graphForecast.writeHTML(outFile);
			
			//print the classic forecast Table
			String tableFile = outFile.getParent() + "/" + "Table.html"; //file must be named Table.html
			try {
				System.out.println("Saving forecast Table to: " + tableFile);
				graphForecast.writeHTMLTable(new File(tableFile));
			} catch (Exception e) {
				System.err.println("Couldn't save forecast Table to: " + tableFile);
			}
			
			
			// print selected figures: Expected number distributions
			for (int i = 0; i < aftershockExpectedNumGraph.size(); i++){
				String file = aftershockExpectedNumGraph.get(i).getName().replaceAll("[^a-zA-Z]",  "").toLowerCase();
				file = outFile.getParent() + "/" + file + ".png"; //file must have this name
				try {
					System.out.println("Saving forecastMFD to: " + file);
					aftershockExpectedNumGraph.get(i).saveAsPNG(file);
				} catch (Exception e) {
					System.err.println("Couldn't save forecastMFD to: " + file);
				}
			}
			
			// cumulative events with time
			try {
				cmlNumGraph.saveAsPNG(outFile.getParent() + "/forecastCmlNum.png");
			} catch (Exception e) {
				System.err.println("Couldn't save forecastCumulativeNumber to: " + outFile.getParent() + "/forecastCmlNum.png");
			}
	
			// ratemap
			if (forecastMapPane != null) {
				for (int i = 0; i < forecastMapPane.getTabCount(); i++){
					GraphWidget graph = (GraphWidget) forecastMapPane.getComponentAt(i);

					String file = graph.getName().replaceAll("[^a-zA-Z]",  "").toLowerCase();
					if (file.contains("rate")){
						file = outFile.getParent() + "/" + file + ".png";
						try {
							System.out.println("Saving forecast Rate Map to: " + file);
							graph.saveAsPNG(file);
						} catch (Exception e) {
							System.err.println("Couldn't save forecast Rate Map to: " + file);
						}
					} else {
//						System.out.println("Skipping " + file);
					}
				}
			} else {
				System.out.println("No Rate map computed.");
			}
						
			// shaking maps
			if (forecastMapPane != null) {
				for (int i = 0; i < forecastMapPane.getTabCount(); i++){
					GraphWidget graph = (GraphWidget) forecastMapPane.getComponentAt(i);

					String file = graph.getName().replaceAll("[^a-zA-Z]",  "").toLowerCase();
					if (file.contains("shaking")){
						file = outFile.getParent() + "/" + file + ".png";
						try {
							System.out.println("Saving forecast Shaking Map to: " + file);
							graph.saveAsPNG(file);
						} catch (Exception e) {
							System.err.println("Couldn't save forecast Shaking Map to: " + file);
						}
					} else {
//						System.out.println("Skipping " + file);
					}
				}
			} else {
				System.out.println("No shaking maps computed.");
			}

			// Write Contours to file
			for(ForecastDuration foreDur : ForecastDuration.values()){
				// find the matching contour set
				int nmatch = 0;
				if(contourList != null) {
				for (ContourModel contours : contourList){
					if (contours.getName().contains(foreDur.toString())){
						nmatch++;
						String name = "contour-" + foreDur.toString();
						File file = new File(outFile.getParent() + "/" + name + ".kml");
						System.out.println("Saving contours to: " + file);
						ETAS_RateModel2D.writeContoursAsKML(contours.getContours(), contours.getName(), contours.getUnits(), file, contours.getCPT());
					}
				}
				if (nmatch > 1) System.out.println("More than one set of contours found for duration: " + foreDur);
				else if (nmatch == 0) System.out.println("No contours found for duration: " + foreDur);
				} else {
					System.out.println("No forecast contours computed.");
				}
				
			}
			
			
			//write a blank logo to the output directory
			// load the data
			String pngFile = "Logo.png";
			InputStream logoIS = GraphicalForecast.class.getResourceAsStream("resources/" + pngFile);
			if (logoIS != null){
			
				File destination = new File(outFile.getParent() + "/" + pngFile);

				try {
					FileUtils.copyInputStreamToFile(logoIS, destination);
				} catch (IOException e) {
				    e.printStackTrace();
					System.err.println("Couldn't copy: " + pngFile + " to file: " + destination);
				}
			} else {
				System.err.println("Couldn't locate file: " + pngFile);
			}
			
		}
	}

	private void setChangeListenerEnabled(boolean enabled){
		changeListenerEnabled = enabled;
	}

	private void link_alpha(){
		if (bParam.getValue() != null && bParam.getValue() < 1d)
			alphaParam.setValue(bParam.getValue()); //link alpha to b if b<1.
		else
			alphaParam.setValue(1d);

//		if(seqSpecModel != null)
//			seqSpecModel.set_alpha(alphaParam.getValue());
		if(bayesianModel != null)
			bayesianModel.set_alpha(alphaParam.getValue());
	}

	private void resetWorkspace() {
		if(verbose) System.out.println("Resetting workspace...");
		setEnableParamsPostFetch(false);
		setEnableParamsPostForecast(false);
		setEnableParamsPostComputeB(false);
		setEnableParamsPostAftershockParams(false);
		
		mainshock = null;
		genericParams = null;
		region = null;
		aftershocks = null;
		aftershockMND = null;
		ogataMND = null;
		genericModel = null;
		bayesianModel = null;
//		seqSpecModel = null;
		faultTrace = null;
		contourList = null;
		
	}
	
	private void doPostFetchCalculations() {
		// update the grid spacing based on the region size
		mapGridSpacingParam.setValue(Math.max(0.01,
				round((region.getMaxLat() - region.getMinLat())*111.111/20,1) ));
		mapGridSpacingParam.getEditor().refreshParamEditor();
		
		// reset the fit constraint sub menu with range around new generic values
		resetFitConstraints(genericParams);
	}
	
	private void doPostFetchPlots() {
		if(verbose) System.out.println("Generating plots...");
		plotAftershockHypocs();
		
		if (aftershockMND == null) 
			aftershockMND = ObsEqkRupListCalc.getMagNumDist(aftershocks, -3.05, 131, 0.1);
		
		plotMFDs(aftershockMND, mcParam.getValue());
		plotMagVsTime();
		plotCumulativeNum();
		plotCatalogText();
		
		setEnableParamsPostFetch(true);
		setEnableParamsPostComputeB(false);
		setEnableParamsPostAftershockParams(false);
		resetFitConstraints(genericParams);
	}
	
	private void refreshTabs(int selectedTab) {
		if (tabbedPane.getTabCount() >= catalog_tab_index && selectedTab >= catalog_tab_index && !tabbedPane.isEnabledAt(catalog_tab_index)){
			tabbedPane.setEnabledAt(catalog_tab_index, true);
			tabbedPane.setForegroundAt(catalog_tab_index, new Color(0,0,0));
		} else 
			Preconditions.checkState(tabbedPane.getTabCount() > catalog_tab_index, "Plots added out of order");

		if (tabbedPane.getTabCount() >= epicenter_tab_index && selectedTab >= epicenter_tab_index && !tabbedPane.isEnabledAt(epicenter_tab_index)){
			tabbedPane.removeTabAt(epicenter_tab_index);
			tabbedPane.insertTab("Epicenters", null, epicenterGraph, "Epicenter Map", epicenter_tab_index);
		} else
			Preconditions.checkState(tabbedPane.getTabCount() > epicenter_tab_index, "Plots added out of order");

		if (tabbedPane.getTabCount() >= mag_num_tab_index && selectedTab >= mag_num_tab_index && !tabbedPane.isEnabledAt(mag_num_tab_index)){
			tabbedPane.removeTabAt(mag_num_tab_index);
			tabbedPane.insertTab("Mag/Num Dist", null, magNumGraph,
					"Aftershock Magnitude vs Number Distribution", mag_num_tab_index);
		} else
			Preconditions.checkState(tabbedPane.getTabCount() > mag_num_tab_index, "Plots added out of order");

		if (tabbedPane.getTabCount() >= mag_time_tab_index && selectedTab >= mag_time_tab_index && !tabbedPane.isEnabledAt(mag_time_tab_index)){
			tabbedPane.removeTabAt(mag_time_tab_index);
			tabbedPane.insertTab("Mag/Time Plot", null, magTimeGraph,
					"Aftershock Magnitude vs Time Plot", mag_time_tab_index);
		} else
			Preconditions.checkState(tabbedPane.getTabCount() > mag_time_tab_index, "Plots added out of order");

		if (tabbedPane.getTabCount() >= cml_num_tab_index && selectedTab >= cml_num_tab_index && !tabbedPane.isEnabledAt(cml_num_tab_index)){
			tabbedPane.removeTabAt(cml_num_tab_index);
			tabbedPane.insertTab("Cumulative Num Plot", null, cmlNumGraph,
					"Cumulative Number Of Aftershocks Plot", cml_num_tab_index);
		} else
			Preconditions.checkState(tabbedPane.getTabCount() > cml_num_tab_index, "Plots added out of order");

		if (tabbedPane.getTabCount() >= pdf_tab_index && selectedTab >= pdf_tab_index && !tabbedPane.isEnabledAt(pdf_tab_index)){
			tabbedPane.removeTabAt(pdf_tab_index);
			tabbedPane.insertTab("Model PDFs", null, pdfGraphsPane,
					"Aftershock Model Prob Dist Funcs", pdf_tab_index);
		} else
			Preconditions.checkState(tabbedPane.getTabCount() > pdf_tab_index, "Plots added out of order");

		if (tabbedPane.getTabCount() >= aftershock_expected_index && selectedTab >= aftershock_expected_index && !tabbedPane.isEnabledAt(aftershock_expected_index)){
			tabbedPane.removeTabAt(aftershock_expected_index);
			tabbedPane.insertTab("Forecast", null, forecastMFDPane,
					"Forecast Plots", aftershock_expected_index);
		}
		else
			Preconditions.checkState(tabbedPane.getTabCount() > aftershock_expected_index, "Plots added out of order");
	
		if (tabbedPane.getTabCount() >= forecast_table_tab_index && selectedTab >= forecast_table_tab_index && !tabbedPane.isEnabledAt(forecast_table_tab_index)){
			tabbedPane.removeTabAt(forecast_table_tab_index);
			tabbedPane.insertTab("Forecast Table", null, forecastTablePane,
					"Forecast Table", forecast_table_tab_index);
		}
		else
			Preconditions.checkState(tabbedPane.getTabCount() > forecast_table_tab_index, "Plots added out of order");
		
		if (tabbedPane.getTabCount() >= forecast_map_tab_index && selectedTab >= forecast_map_tab_index && !tabbedPane.isEnabledAt(forecast_map_tab_index)){
			tabbedPane.removeTabAt(forecast_map_tab_index);
			tabbedPane.insertTab("Forecast Maps", null, forecastMapPane,
					"Forcast Maps (rate, intensity)", forecast_map_tab_index);
		}
		else
			Preconditions.checkState(tabbedPane.getTabCount() > forecast_map_tab_index, "Plots added out of order");
		
		
		tabbedPane.setSelectedIndex(selectedTab);
	}
	
	private void refreshTimeWindowEditor(){
		boolean now = nowBoolean.getValue();
		
		forecastStartTimeParam.getEditor().setVisible(!now);
		
		if (forecastStartTimeParam.getValue() == null){
			fetchButton.getEditor().setEnabled(now);
			fetchButton.getEditor().refreshParamEditor();
			quickForecastButton.getEditor().setEnabled(now);
			quickForecastButton.getEditor().refreshParamEditor();
		}
		else{
			fetchButton.getEditor().setEnabled(true);
			fetchButton.getEditor().refreshParamEditor();
			quickForecastButton.getEditor().setEnabled(true);
			quickForecastButton.getEditor().refreshParamEditor();
		}
		
		dataStartTimeParam.getEditor().setEnabled(!now);
		dataEndTimeParam.getEditor().setEnabled(!now);
		forecastStartTimeParam.getEditor().setEnabled(!now);
	}
		
	
	private void resetFitConstraints(GenericETAS_Parameters params){
		if(verbose)System.out.println("Updating fit constraints based on new generic params...");

		boolean restrictParameters = false; 
		
		double new_ams = params.get_ams();
		double new_a = params.get_a();
		double new_p = params.get_p();
		double new_c = params.get_c();
			
		//correct the a-value if Mc is not the same as refMag
		double refMag = params.get_refMag();
		double maxMag = params.get_maxMag();
		double mc = mcParam.getValue();
			
		new_ams += Math.log10((maxMag - refMag)/(maxMag - mc));
		new_a += Math.log10((maxMag - refMag)/(maxMag - mc));
		
		double min_ams = amsValRangeParam.getValue().getLowerBound();
		double max_ams = amsValRangeParam.getValue().getUpperBound();
		
		double min_a = aValRangeParam.getValue().getLowerBound();
		double max_a = aValRangeParam.getValue().getUpperBound();
		
		double min_p = pValRangeParam.getValue().getLowerBound();
		if (min_p < 0.25)
			min_p = 0.25;
				
		double max_p = pValRangeParam.getValue().getUpperBound();
		if (max_p > 2)
			max_p = 2;
		
		double min_c = cValRangeParam.getValue().getLowerBound();
		if (min_c <= 1e-5)
			min_c = 1e-5;
		
		double max_c = cValRangeParam.getValue().getUpperBound();
		if(max_c > 1)
			max_c = 1;
		if(max_c < min_c)
			max_c = min_c;
		
		// reset the b-value and Mc //no, do this in set magComplete()
//		bParam.setValue(params.get_b());
//		bParam.getEditor().refreshParamEditor();
//		mcParam.setValue(mc);
//		mcParam.getEditor().refreshParamEditor();
		
		// recenter the ranges around the new parameter values
		amsValRangeParam.setValue(new Range(round(new_ams - (max_ams - min_ams)/2,2), round(new_ams + (max_ams - min_ams)/2,2)));
		amsValRangeParam.getEditor().refreshParamEditor();
		
		if(restrictParameters && !devMode && aftershocks.getRupsAboveMag(mc).size() < 10){
			aValRangeParam.setValue(new Range(round(new_a,2), round(new_a,2)));
			aValRangeParam.getEditor().refreshParamEditor();
			aValNumParam.setValue(1);
			aValNumParam.getEditor().refreshParamEditor();
			aValRangeParam.getEditor().setEnabled(false);
			aValNumParam.getEditor().setEnabled(false);
		} else {	
			aValRangeParam.setValue(new Range(round(new_a - (max_a - min_a)/2,2), round(new_a + (max_a - min_a)/2,2)));
			aValRangeParam.getEditor().refreshParamEditor();
			aValNumParam.getEditor().refreshParamEditor();
			aValRangeParam.getEditor().setEnabled(true);
			aValNumParam.getEditor().setEnabled(true);
		}
		
		if(restrictParameters && !devMode && aftershocks.getRupsAboveMag(mc).size() < 20){
			pValRangeParam.setValue(new Range(round(new_p,2), round(new_p,2)));
			pValRangeParam.getEditor().refreshParamEditor();
			pValNumParam.setValue(1);
			pValNumParam.getEditor().refreshParamEditor();
			pValNumParam.getEditor().setEnabled(false);
			pValRangeParam.getEditor().setEnabled(false);
		} else {
			pValRangeParam.setValue(new Range(round(new_p - (max_p - min_p)/2,2), round(new_p + (max_p - min_p)/2,2)));
			pValRangeParam.getEditor().refreshParamEditor();
			pValNumParam.getEditor().refreshParamEditor();
			pValNumParam.getEditor().setEnabled(true);
			pValRangeParam.getEditor().setEnabled(true);
		}

		if(restrictParameters && !devMode && aftershocks.getRupsAboveMag(mc).size() < 30){
			cValRangeParam.setValue(new Range(round(new_c,6), round(new_c,6)));
			cValRangeParam.getEditor().refreshParamEditor();
			cValNumParam.setValue(1);
			cValNumParam.getEditor().refreshParamEditor();
			cValNumParam.getEditor().setEnabled(false);
			cValRangeParam.getEditor().setEnabled(false);
		} else {
			cValRangeParam.setValue(new Range(round(new_c / Math.sqrt(max_c/min_c),6), round(new_c * Math.sqrt(max_c/min_c),6)));
			cValRangeParam.getEditor().refreshParamEditor();
			cValNumParam.getEditor().refreshParamEditor();
			cValNumParam.getEditor().setEnabled(true);
			cValRangeParam.getEditor().setEnabled(true);
		}

		
		// checks ranges and updates the num parameter based on the range
		updateRangeParams(amsValRangeParam, amsValNumParam, 51);
		updateRangeParams(aValRangeParam, aValNumParam, 51);
		updateRangeParams(pValRangeParam, pValNumParam, 51);
		updateRangeParams(cValRangeParam, cValNumParam, 51);
		
		//clear the sequence specific and Bayesian fits
//		if (seqSpecModel != null)
//			seqSpecModel = null;
		if (bayesianModel != null)
			bayesianModel = null;
		
		
		amsValParam.setValue(null);
		amsValParam.getEditor().refreshParamEditor();
		aValParam.setValue(null);
		aValParam.getEditor().refreshParamEditor();
		pValParam.setValue(null);
		pValParam.getEditor().refreshParamEditor();
		cValParam.setValue(null);
		cValParam.getEditor().refreshParamEditor();
		
		tabbedPane.setForegroundAt(pdf_tab_index, new Color(128,128,128));
		tabbedPane.setEnabledAt(pdf_tab_index, false);
		
		tabbedPane.setForegroundAt(aftershock_expected_index, new Color(128,128,128));
		tabbedPane.setEnabledAt(aftershock_expected_index, false);

		tabbedPane.setForegroundAt(forecast_table_tab_index, new Color(128,128,128));
		tabbedPane.setEnabledAt(forecast_table_tab_index, false);

		tabbedPane.setForegroundAt(forecast_map_tab_index, new Color(128,128,128));
		tabbedPane.setEnabledAt(forecast_map_tab_index, false);
		

//		tabbedPane.setSelectedIndex(mag_num_tab_index); //make sure to look away first!
		
		// need to remove forecast from cumulativeNumber plot
	}
	
	private void updateRangeParams(RangeParameter rangeParam, IntegerParameter numParam, int defaultNum) {
		Preconditions.checkState(defaultNum > 1);
		Range range = rangeParam.getValue();
		if (range == null)
			return;
		boolean same = range.getLowerBound() == range.getUpperBound();
		if (same && numParam.getValue() > 1)
			numParam.setValue(1);
		else if (!same && numParam.getValue() == 1)
			numParam.setValue(defaultNum);
		numParam.getEditor().refreshParamEditor();
	}
	
	/**
	 * disables/enables all parameters that are dependent on the fetch step and beyond
	 */
	private void setEnableParamsPostFetch(boolean enabled) {
	
		saveCatalogButton.getEditor().setEnabled(enabled);
		computeBButton.getEditor().setEnabled(enabled);
		tectonicRegimeParam.getEditor().setEnabled(enabled);
		
		// these used to be enabled after computing b but we now allow the user to just use default B
		
		mcParam.getEditor().setEnabled(enabled);
		if (!enabled || (enabled && devMode)){
			bParam.getEditor().setEnabled(enabled);
			magPrecisionParam.getEditor().setEnabled(enabled);
		}
		
		amsValRangeParam.getEditor().setEnabled(enabled);
		amsValNumParam.getEditor().setEnabled(enabled);
		aValRangeParam.getEditor().setEnabled(enabled);
		aValNumParam.getEditor().setEnabled(enabled);
		pValRangeParam.getEditor().setEnabled(enabled);
		pValNumParam.getEditor().setEnabled(enabled);
		cValRangeParam.getEditor().setEnabled(enabled);
		cValNumParam.getEditor().setEnabled(enabled);
		computeAftershockParamsButton.getEditor().setEnabled(enabled);
		computeAftershockForecastButton.getEditor().setEnabled(enabled);
		generateMapButton.getEditor().setEnabled(enabled);
		
		rmaxParam.getEditor().setEnabled(enabled && timeDepMcParam.getValue());
	}
	
	/**
	 * disables all parameters that are dependent on the compute b step and beyond
	 */
	private void setEnableParamsPostComputeB(boolean enabled) {
		if (!enabled){
			setEnableParamsPostAftershockParams(enabled);

			// disable tabs
			for (int i = aftershock_expected_index; i < tabbedPane.getTabCount() ; i++){
				tabbedPane.setForegroundAt(i, new Color(128,128,128));
				tabbedPane.setEnabledAt(i, false);
			}
		}
	}

	private void setEnableParamsPostComputeMc(boolean enabled) {
		if (!enabled){
			setEnableParamsPostAftershockParams(enabled);

			for (int i = pdf_tab_index; i < tabbedPane.getTabCount() ; i++){
				tabbedPane.setForegroundAt(i, new Color(128,128,128));
				tabbedPane.setEnabledAt(i, false);
			}
		}
	}
	
	private void setEnableParameterEditing(boolean enabled){ 
		//	 these should not be editable in safeMode, but should be responsive in devMode	
		amsValParam.getEditor().setEnabled(enabled); // no capability to set directly yet (add in, for custom forecast?)
		aValParam.getEditor().setEnabled(enabled); // no capability to set directly yet
		pValParam.getEditor().setEnabled(enabled); // no capability to set directly yet
		cValParam.getEditor().setEnabled(enabled); // no capability to set directly yet
//		if (!enabled) {
//			amsValParam.getEditor().getComponent().setForeground(Color.gray);
//			aValParam.getEditor().getComponent().setForeground(Color.lightGray);
//			pValParam.getEditor().getComponent().setForeground(Color.lightGray);
//			cValParam.getEditor().getComponent().setForeground(Color.lightGray);
//		}
			
	}

	private void setEnableParamsPostAftershockParams(boolean enabled) {
		
		if (enabled) {
			computeAftershockForecastButton.setButtonText("Run Specific Forecast");
			computeAftershockForecastButton.setInfo("Compute aftershock forecast using parameters estimated for this sequence");
		} else {
			computeAftershockForecastButton.setButtonText("Run Generic Forecast");
			computeAftershockForecastButton.setInfo("Compute aftershock forecast using typical parameters");
		}
		computeAftershockForecastButton.getEditor().refreshParamEditor();
		
		for (int i = aftershock_expected_index; i < tabbedPane.getTabCount() ; i++){
			tabbedPane.setForegroundAt(i, new Color(128,128,128));
			tabbedPane.setEnabledAt(i, false);
		}
		
		plotSpecificOnlyParam.getEditor().setEnabled(enabled);
		if (!enabled)
			plotSpecificOnlyParam.setValue(false);
		plotSpecificOnlyParam.getEditor().refreshParamEditor();
		
		setEnableParamsPostForecast(false);
//		if (!enabled){
//			seqSpecModel = null;
//		}
	}
	
	private void setEnableParamsPostForecast(boolean enabled) {
		if (enabled)
			for (int i = forecast_map_tab_index; i < tabbedPane.getTabCount() ; i++){
				tabbedPane.setForegroundAt(i, new Color(128,128,128));
				tabbedPane.setEnabledAt(i, false);
			}
		
		generateMapButton.getEditor().setEnabled(enabled);
		generateMapButton.getEditor().refreshParamEditor();
		
		setEnableParamsPostRender(false); //either way, reset the plot
	}
	
	private void setEnableParamsPostRender(boolean enabled) {
		publishAdvisoryButton.getEditor().setEnabled(enabled);
		publishAdvisoryButton.getEditor().refreshParamEditor();
		
		if (!enabled) {
			//reset the map curves if any
			pgvCurves=null;
			pgaCurves=null;
			psaCurves=null;
		}
	}
	
	private void validateRange(Range range, int num, String name) {
		Preconditions.checkState(range != null, "Must supply "+name+" range");
		boolean same = range.getLowerBound() == range.getUpperBound();
		if (same)
			Preconditions.checkState(num == 1, "Num must equal 1 for fixed "+name);
		else
			Preconditions.checkState(num > 1, "Num must be >1 for variable "+name);
	}

	private static void validateParameter(Double value, String name) {
		Preconditions.checkState(value != null, "Must specify "+name);
		Preconditions.checkState(Doubles.isFinite(value), name+" must be finite: %s", value);
	}

	private void checkArguments(String... args){
		//set defaults
		verbose = false;
		devMode = false;
		D = false;
			
		//check for arguments
		for (String argument : args){
			if (argument.contains("verbose")) verbose = true;
			if (argument.contains("expert")) devMode = true;
			if (argument.contains("debug")) D = true;
			if(verbose) System.out.println(argument);
		}
	}
	
	private Location getCentroid() {
		return ETAS_StatsCalc.getCentroid(mainshock, aftershocks);
	}

	private double getTimeSinceMainshock(ObsEqkRupture rup) {
		long ms = mainshock.getOriginTime();
		long as = rup.getOriginTime();
		long delta = as - ms;
		return (double)delta/(1000*60*60*24);
	}

	private double getTimeRemainingInUTCDay(){
		double daysLeftInDay;
		TimeZone.setDefault(utc);
		GregorianCalendar origin = mainshock.getOriginTimeCal();
		
		SimpleDateFormat formatter=new SimpleDateFormat("d MMM yyyy, HH:mm:ss");
		formatter.setTimeZone(utc); //utc=TimeZone.getTimeZone("UTC"));
		
		GregorianCalendar daybreak = new GregorianCalendar(
				origin.get(GregorianCalendar.YEAR), origin.get(GregorianCalendar.MONTH), origin.get(GregorianCalendar.DAY_OF_MONTH));
//		daybreak.setTimeZone(origin.getTimeZone());
		
		if(D) System.out.println(formatter.format(origin.getTime()));
		if(D) System.out.println(formatter.format(daybreak.getTime()));
		
		daysLeftInDay = 1 - (double) (origin.getTimeInMillis() - daybreak.getTimeInMillis())/ETAS_StatsCalc.MILLISEC_PER_DAY;
		return daysLeftInDay;
	}
//	
	
	
	private double round(double val){
		return round(val, sigDigits);
		
	}
	private double round(double val, int sigDigits){
		return Math.round(val*Math.pow(10, sigDigits))/Math.pow(10, sigDigits);
	}

	private void loadSourceFromFile() {
//		File finiteSourceFile = new File();

		String source = " 34.6985 -118.5089\n"
				+"34.5478 -118.1039\n"
				+"34.4029 -117.7536\n"
				+"34.3163 -117.5490\n"
				+"34.3163 -117.5490\n"
				+"34.2709 -117.4510\n"
				+"34.2328 -117.3887\n"
				+"34.1731 -117.2742\n"
				+"34.1500 -117.2220\n"
				+"34.1500 -117.2220\n"
				+"34.0928 -117.0677\n"
				+"34.0738 -117.0139\n"
				+"34.0338 -116.9023\n"
				+"34.0113 -116.8735\n"
				+"33.9591 -116.8198\n"
				+"33.9532 -116.8014\n"
				+"33.9374 -116.7786\n"
				+"33.9442 -116.6858\n"
				+"33.9176 -116.6239\n"
				+"33.9070 -116.5849\n"
				+"33.8847 -116.5169\n"
				+"33.8481 -116.4265\n"
				+"33.8485 -116.3830\n"
				+"33.7882 -116.2463\n"
				+"33.7882 -116.2463\n"
				+"33.3501 -115.7119";
		
//		List<String> lines = source.split("\\n");
		String[] lines = source.split("\\n");
		
		if (faultTrace == null)
			faultTrace = new FaultTrace("custom");
	
		for (int i=0; i<lines.length; i++) {
			String line = lines[i];
			if (line.startsWith("#")) {
				//header
			} else {
				if (!line.isEmpty()) {
					String[] pts = line.trim().split(" ");
					if(D) System.out.println(i + ": " + line +  " " + pts[0] + " " + pts[1]);
					double lat = Double.parseDouble(pts[0]);
					double lon = Double.parseDouble(pts[1]);

					
					faultTrace.add(new Location(lat, lon));
				}
			}
		}
		System.out.println("Loaded "+faultTrace.size()+" points from finite source file.");
	}

	private GeoFeatureList loadCities(Region region, int number){
		Boolean D = false;
		
		GeoFeatureList cities = new GeoFeatureList();
		if(D) System.out.println("Loading Geographic information...");

		// load the data
		InputStream citiesIS = GeoFeatureList.class.getResourceAsStream("resources/worldcities1000.txt");
		List<String> lines = new ArrayList<String>();
		try{
			lines = IOUtils.readLines(citiesIS, StandardCharsets.UTF_8);
		} catch (IOException e) {
			System.out.println("Couldn't load city information");
		}

		//populate the feature list
		for (String line: lines){
			cities.addFeatureFromLine(line);
		}
		
		if(D) System.out.println(cities.size() + " cities added to list.");
		
		cities.getFeaturesInRegion(region);
//		cities.getFeaturesAboveLevel(minLevel);
		cities.sortByMapLevel(false); //false --> descending order
		cities.thinFeatures(Math.sqrt(region.getExtent())/10);

		if (cities.size() > number)
			cities.getFeaturesAboveLevel(cities.get(number).mapLevel);
		

		if(D)
			for(GeoFeature city : cities)
				System.out.println(city);
		
		if(D) System.out.println(cities.size() + " cities in mapped region.");
		
		return cities;
	}
	
	private ArrayList<XY_DataSet> loadCountries(Region region){
		if(verbose) System.out.println("Loading country borders...");
		
		// Load country border lat/long points
		InputStream bordersIS = GeoFeatureList.class.getResourceAsStream("resources/CountryBorders.txt");
			
		List<String> lines = new ArrayList<String>();
		try{
			lines = IOUtils.readLines(bordersIS, StandardCharsets.UTF_8);
		} catch (IOException e) {
			System.out.println("Couldn't load country border information");
		}

		ArrayList<XY_DataSet> countryBorders = new ArrayList<XY_DataSet>();
//		ArrayList<String> countryNames = new ArrayList<String>();
		for (int i=1; i<lines.size(); i++) { // Skip first line of file (headers)
			String line = lines.get(i);
			String[] elem = line.split(", ");
			String name = elem[0];
			int numberPoints =  Integer.parseInt(elem[1]);
			XY_DataSet points = new DefaultXY_DataSet();
			
			Boolean isInside = false;
			for (int j=0; j<numberPoints; j++) {
				float lat = Float.parseFloat(elem[2*j+2]);
				float lon = Float.parseFloat(elem[2*j+3]);
				Location loc = new Location(lat,lon);
				points.set(lon, lat);
				if (region.contains(loc))
					isInside = true;
			}	
		
			if(isInside) {
				points.setName(name);
				countryBorders.add(points); 
			}
		}
		
		
		return countryBorders;
	}
	
	private GriddedGeoDataSet cloneGriddedGeoDataSet(GriddedGeoDataSet dataSet) {
		GriddedGeoDataSet newDataSet = new GriddedGeoDataSet(dataSet.getRegion(), dataSet.isLatitudeX());
		for (int i = 0; i < dataSet.size(); i++) {
			newDataSet.set(i, dataSet.get(i));
		}
		return newDataSet;
	}
	
	
	private List<String> tipText;
	private void printTip(int step){
		
		if(tipText == null){
			tipText = new ArrayList<String>();
		}
		
		SimpleDateFormat formatter=new SimpleDateFormat("d MMM yyyy");
		formatter.setTimeZone(utc); //utc=TimeZone.getTimeZone("UTC"));
		double elapsedDays = (double) (System.currentTimeMillis() - expirationDate.getTimeInMillis())/ETAS_StatsCalc.MILLISEC_PER_DAY;
		
		String welcomeMessage = "This a Beta version of the Aftershock Forecaster software. Get the latest version from www.opensha.org/apps.\n"
				+ "The Beta version will expire " + formatter.format(expirationDate.getTime()) + String.format(" (%d days remaining).", (int) -elapsedDays);

		tipText.add(welcomeMessage + "\n\n>> Welcome to the aftershock forecaster. Enter a USGS event ID to get started. (e.g: ci14607652)");
		tipText.add(">> Specify a forecast start time and duration. Then click \"Fetch Data\" to retrieve the catalog\n  ...or click \"Quick Forecast\" to run the entire forecast automatically with default settings.");
		tipText.add(">> Click \"Compute Model Fit\" to compute sequence-specific model\n  ...or go straight to \"Run Generic Forecast\" to get a generic forecast for this region.");
		tipText.add(">> Click \"Compute Model Fit\" to compute sequence-specific model\n  ...or click \"Render\" to get a generic aftershock rate map.");
		tipText.add(">> Click \"Run Specific Forecast\" to get a sequence-specific forecast.");
		tipText.add(">> Choose the map settings and click \"Render\" to generate forecast maps.");
		tipText.add(">> Go back and edit the forecast or click \"Publish\" to generate a printable forecast document.");
		tipText.add(">> Click \"Run Specific Forecast\" to (re)generate forecast plots.");
		tipText.add(">> Click \"Render\" to (re)generate the spatial forecast plots.");
		if (tipsOn) System.out.println(tipText.get(step));
	}
	
	public static void main(String... args) {
		
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run(){
				new AftershockStatsGUI_ETAS(args);
			}
		});

		
	}

}
