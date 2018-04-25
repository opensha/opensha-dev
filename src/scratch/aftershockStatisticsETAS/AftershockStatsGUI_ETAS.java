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
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.EnumSet;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.SimpleTimeZone;
import java.util.TimeZone;
import java.util.concurrent.TimeUnit;

import javax.swing.*;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.dom4j.Element;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYDrawableAnnotation;
import org.jfree.chart.labels.XYItemLabelGenerator;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.data.Range;
import org.jfree.data.xy.AbstractXYDataset;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.RectangleEdge;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.function.AbstractXY_DataSet;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DatasetBinner;
import org.opensha.commons.data.siteData.impl.TectonicRegime;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GeoTools;
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
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.param.constraint.impl.DoubleConstraint;
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
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupListCalc;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.CalcProgressBar;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

import wContour.Contour;
import wContour.Global.PointD;
import wContour.Global.PolyLine;

import java.awt.Font;
import java.awt.Image;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;


public class AftershockStatsGUI_ETAS extends JFrame implements ParameterChangeListener {
		
	public AftershockStatsGUI_ETAS(String... args) {
		checkArguments(args);
		if (D) System.out.println("verbose = " + verbose + ", expert mode = " + expertMode);
		createAndShowGUI();
	}
    
	private boolean D = true; //debug
	private boolean expertMode;
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
	
	private IntegerParameter plotMMIParam;
	private DoubleParameter gridSpacingParam;
	private DoubleParameter plotPOEParam;
	
	private DoubleParameter forecastStartTimeParam;
	private DoubleParameter forecastEndTimeParam;
	private DoubleParameter forecastDurationParam;
	
	private ButtonParameter computeAftershockForecastButton;
	private ButtonParameter generateMapButton;
	private BooleanParameter fitMainshockSourceParam;
	private BooleanParameter fitShakeMapSourceParam;
	
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
//	private GraphWidget aftershockExpectedNumGraph;
//	private GraphWidget aftershockProbabilityGraph;
	
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
	
	private ComcatAccessor accessor;
	private WC1994_MagLengthRelationship wcMagLen;
	
	private Region region;
	private ObsEqkRupture mainshock;
	private FaultTrace faultTrace;
	private ObsEqkRupList aftershocks;
	private ObsEqkRupture largestShock;
	private ETAS_RateModel2D rateModel2D;
	private List<ContourModel> contourList;
	private String shakeMapURL;
	
	private OgataMagFreqDist ogataMND;
		
	private GenericETAS_ParametersFetch genericFetch = null;
	private GenericETAS_Parameters genericParams = null;
	private ETAS_AftershockModel_Generic genericModel = null;
	private ETAS_AftershockModel_SequenceSpecific seqSpecModel = null;
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
//	private int consoleHeight = (int) (80); 
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
		
		eventIDParam = new StringParameter("USGS Event ID");
		eventIDParam.setValue("us20002926");
		eventIDParam.setInfo("Get IDs from https://earthquake.usgs.gov/earthquakes/");
		eventIDParam.addParameterChangeListener(this);
		forecastParams.addParameter(eventIDParam);
		
		nowBoolean = new BooleanParameter("Begin at current time", true);
		nowBoolean.setInfo("Use current time for forecast start");
		nowBoolean.addParameterChangeListener(this);
		forecastParams.addParameter(nowBoolean);
		
		forecastStartTimeParam = new DoubleParameter("Forecast Start Time", 0d, 3660);
		forecastStartTimeParam.setUnits("Days");
		forecastStartTimeParam.addParameterChangeListener(this);
		forecastParams.addParameter(forecastStartTimeParam);

		forecastDurationParam = new DoubleParameter("Forecast Duration", 0d, 366, new Double(366d));
		forecastDurationParam.setUnits("Days");
		forecastDurationParam.addParameterChangeListener(this);
		forecastParams.addParameter(forecastDurationParam);

		plotAllDurationsParam = new BooleanParameter("Plot Week/Month/Year", true);	
		plotAllDurationsParam.addParameterChangeListener(this);
		forecastParams.addParameter(plotAllDurationsParam);
		
		quickForecastButton = new ButtonParameter("Forecast using default settings", "Quick Forecast");
		quickForecastButton.setInfo("Forecast using default settings");
		quickForecastButton.addParameterChangeListener(this);
		quickForecastButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		forecastParams.addParameter(quickForecastButton);
		
		timeWindow = new ParameterList(); 
		dataStartTimeParam = new DoubleParameter("Data Start Time", 0d, 366, new Double(0d));
		dataStartTimeParam.setUnits("Days");
		dataStartTimeParam.setInfo("Relative to main shock origin time");
		dataStartTimeParam.addParameterChangeListener(this);
		timeWindow.addParameter(dataStartTimeParam);
		
		dataEndTimeParam = new DoubleParameter("Data End Time", 0d, 3660);
		dataEndTimeParam.setUnits("Days");
		dataEndTimeParam.setInfo("Relative to main shock origin time");
		dataEndTimeParam.addParameterChangeListener(this);
		timeWindow.addParameter(dataEndTimeParam);
				
		forecastEndTimeParam = new DoubleParameter("Forecast End Time", 0d, 3660+366);
		forecastEndTimeParam.setUnits("Days");
//		forecastEndTimeParam.addParameterChangeListener(this);
//		timeWindow.addParameter(forecastEndTimeParam);
		
		timeWindowEditParam = new ParameterListParameter("Edit data time window", timeWindow);
		if (expertMode) forecastParams.addParameter(timeWindowEditParam);		
		
		// these are inside region editor
		regionTypeParam = new EnumParameter<AftershockStatsGUI_ETAS.RegionType>(
				"Aftershock Zone Type", EnumSet.allOf(RegionType.class), RegionType.CIRCULAR_WC94, null);
		regionTypeParam.setInfo("For collecting aftershocks");
		regionTypeParam.addParameterChangeListener(this);
//		dataParams.addParameter(regionTypeParam);

		regionCenterLatParam = new DoubleParameter("Center Latitude",
				new DoubleConstraint(GeoTools.LAT_MIN,GeoTools.LAT_MAX),
				DECIMAL_DEGREES, null);
		regionCenterLatParam.getConstraint().setNullAllowed(true);
		regionCenterLonParam = new DoubleParameter("Center Longitude",
				new DoubleConstraint(GeoTools.LON_MIN,360),
				DECIMAL_DEGREES, null);
		regionCenterLonParam.getConstraint().setNullAllowed(true);
		
		radiusParam = new DoubleParameter("Radius", 0d, 1000, new Double(20));
		radiusParam.setUnits("km");
		minLatParam = new DoubleParameter("Min Lat", -90d, 90d, new Double(0));
		maxLatParam = new DoubleParameter("Max Lat", -90d, 90d, new Double(0));
		minLonParam = new DoubleParameter("Min Lon", -180d, 360d, new Double(0));
		maxLonParam = new DoubleParameter("Max Lon", -180d, 360d, new Double(0));
		
		minDepthParam = new DoubleParameter("Min Depth", 0d, 1000d, new Double(0));
		minDepthParam.setUnits("km");
		maxDepthParam = new DoubleParameter("Max Depth", 0d, 1000d, new Double(1000d));
		maxDepthParam.setUnits("km");
		
		regionCenterTypeParam = new EnumParameter<AftershockStatsGUI_ETAS.RegionCenterType>(
				"Aftershock Zone Center", EnumSet.allOf(RegionCenterType.class), RegionCenterType.CENTROID, null);

//		regionCenterLocParam = new LocationParameter("Aftershock Zone Center Location");

		regionCenterTypeParam.addParameterChangeListener(this);
		
		regionList = new ParameterList();
		
		regionEditParam = new ParameterListParameter("Edit Aftershock Zone", regionList);
		regionEditParam.setInfo("To set location constraints");
		
		regionEditParam.addParameterChangeListener(this);
		regionEditParam.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		dataParams.addParameter(regionEditParam);
		
		fetchButton = new ButtonParameter("USGS Event Webservice", "Fetch Data");
		fetchButton.setInfo("From USGS ComCat");
		fetchButton.addParameterChangeListener(this);
		fetchButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		dataParams.addParameter(fetchButton);
		
		loadCatalogButton = new ButtonParameter("External Catalog", "Load Local Catalog");
		loadCatalogButton.setInfo("Load catalog in 10 column format");
		loadCatalogButton.addParameterChangeListener(this);
		loadCatalogButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		dataParams.addParameter(loadCatalogButton);
		
		saveCatalogButton = new ButtonParameter("Aftershock Catalog", "Save Local Catalog");
		saveCatalogButton.setInfo("Save catalog in 10 column format");
		saveCatalogButton.addParameterChangeListener(this);
		saveCatalogButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		dataParams.addParameter(saveCatalogButton);
		
		
		/*
		 * Constraint params
		 */
		computeBButton = new ButtonParameter("b-value", "Compute b");
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
		
		if (!expertMode){
			bParam.getEditor().setEnabled(expertMode);
//			mcParam.getEditor().setEnabled(expertMode);
			magPrecisionParam.getEditor().setEnabled(expertMode);
		}
		
		outputParams.addParameter(bParam);
		outputParams.addParameter(mcParam);
		outputParams.addParameter(magPrecisionParam);
		
		
		
//		mfdParams.addParameter(bParam);
//		mfdParams.addParameter(mcParam);
//		mfdParams.addParameter(magPrecisionParam);
		
		timeDepMcParam = new BooleanParameter("Apply time dep. Mc", false);
		timeDepMcParam.setInfo("Apply time dependent magnitude of completeness");
		if (expertMode){
			timeDepMcParam.addParameterChangeListener(this);
			fitParams.addParameter(timeDepMcParam);
		}
		
		alphaParam = new DoubleParameter("alpha-value", 1d);
		alphaParam.setInfo("Linked to b-value");
		alphaParam.addParameterChangeListener(this);
//		dataParams.addParameter(alphaParam);
		
		/* make this display more useful labels */
		tectonicRegimeParam = new EnumParameter<TectonicRegime>(
				"Tectonic Regime", EnumSet.allOf(TectonicRegime.class), TectonicRegime.GLOBAL_AVERAGE, null);
		tectonicRegimeParam.addParameterChangeListener(this);
		fitParams.addParameter(tectonicRegimeParam);
		
		
		//these are inside constraint editor
		fitMSProductivityParam = new BooleanParameter("Fit MS Productivity", true);
		fitMSProductivityParam.addParameterChangeListener(this);
		
		amsValRangeParam = new RangeParameter("ams-value range", new Range(-4.0, -1));
		amsValRangeParam.addParameterChangeListener(this);
		
		amsValNumParam = new IntegerParameter("ams-value num", 1, 101, new Integer(51));
		amsValNumParam.addParameterChangeListener(this);
		
		aValRangeParam = new RangeParameter("a-value range", new Range(-4.0, -1));
		aValRangeParam.addParameterChangeListener(this);
		
		aValNumParam = new IntegerParameter("a-value num", 1, 101, new Integer(31));
		aValNumParam.addParameterChangeListener(this);
		
		pValRangeParam = new RangeParameter("p-value range", new Range(0.5, 2.0));
		pValRangeParam.addParameterChangeListener(this);
		
		pValNumParam = new IntegerParameter("p-value num", 1, 101, new Integer(31));
		pValNumParam.addParameterChangeListener(this);
		
		cValRangeParam = new RangeParameter("c-value range", new Range(1e-5, 1));
		cValRangeParam.addParameterChangeListener(this);
		
		cValNumParam = new IntegerParameter("c-value num", 1, 101, new Integer(21));
		cValNumParam.addParameterChangeListener(this);
		
		rmaxParam = new DoubleParameter("rmax", 1d, Double.POSITIVE_INFINITY, new Double(200));
		rmaxParam.addParameterChangeListener(this);
		
		constraintList = new ParameterList();
		if(expertMode) constraintList.addParameter(fitMSProductivityParam);
		constraintList.addParameter(amsValRangeParam);
		constraintList.addParameter(amsValNumParam);
		constraintList.addParameter(aValRangeParam);
		constraintList.addParameter(aValNumParam);
		constraintList.addParameter(pValRangeParam);
		constraintList.addParameter(pValNumParam);
		constraintList.addParameter(cValRangeParam);
		constraintList.addParameter(cValNumParam);
		if(expertMode) constraintList.addParameter(rmaxParam);
				
		constraintEditParam = new ParameterListParameter("Edit Fit Constraints", constraintList);
		constraintEditParam.setInfo("To set fit constraints");
		constraintEditParam.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		constraintEditParam.addParameterChangeListener(this);
		fitParams.addParameter(constraintEditParam);
				
		computeAftershockParamsButton = new ButtonParameter("Aftershock Params", "Compute Model Fit");
		computeAftershockParamsButton.addParameterChangeListener(this);
		computeAftershockParamsButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		fitParams.addParameter(computeAftershockParamsButton);
		
		/*
		 * Fit params
		 */
		amsValParam = new DoubleParameter("ams", new Double(0d));
		amsValParam.setValue(null);
		if (!expertMode) amsValParam.setNonEditable();
		amsValParam.addParameterChangeListener(this);
		amsValParam.getEditor().getComponent().setMinimumSize(null);
		amsValParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		outputParams.addParameter(amsValParam);
		
		aValParam = new DoubleParameter("a", new Double(0d));
		aValParam.setValue(null);
		aValParam.addParameterChangeListener(this);
		aValParam.getEditor().getComponent().setMinimumSize(null);
		aValParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		outputParams.addParameter(aValParam);
		
		pValParam = new DoubleParameter("p", new Double(0d));
		pValParam.setValue(null);
		if (!expertMode) pValParam.setNonEditable();
		pValParam.addParameterChangeListener(this);
		pValParam.getEditor().getComponent().setMinimumSize(null);
		pValParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		outputParams.addParameter(pValParam);
		
		cValParam = new DoubleParameter("c", new Double(0d));
		cValParam.setValue(null);
		if (!expertMode) cValParam.setNonEditable();
		cValParam.addParameterChangeListener(this);
		cValParam.getEditor().getComponent().setMinimumSize(null);
		cValParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		outputParams.addParameter(cValParam);
		
		plotMMIParam = new IntegerParameter("MMI", 2, 10, new Integer(6) );
		plotMMIParam.setInfo("MMI level for map");
		plotMMIParam.addParameterChangeListener(this);
		plotMMIParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		mapPlotParams.addParameter(plotMMIParam);
		
		plotPOEParam = new DoubleParameter("POE", 0, 1, new Double(0.5));
		plotPOEParam.setInfo("Probability of exceedence level for map");
		plotPOEParam.addParameterChangeListener(this);
		plotPOEParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		mapPlotParams.addParameter(plotPOEParam);
				
		gridSpacingParam = new DoubleParameter("\u0394 (km)", 1, 100, new Double(10));
		gridSpacingParam.setInfo("Cell size for map (km)");
		gridSpacingParam.addParameterChangeListener(this);
		gridSpacingParam.getEditor().getComponent().setPreferredSize(new Dimension(outputWidth, 50));
		mapPlotParams.addParameter(gridSpacingParam);
		
		magRefParam = new DoubleParameter("Reference Magnitude", new Double(0d));
		
		computeAftershockForecastButton = new ButtonParameter("Aftershock Forecast", "Run Generic Forecast");
		computeAftershockForecastButton.addParameterChangeListener(this);
		computeAftershockForecastButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		fitParams.addParameter(computeAftershockForecastButton);
		
		generateMapButton = new ButtonParameter("Forecast Map", "Render");
		generateMapButton.addParameterChangeListener(this);
		generateMapButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		mapParams.addParameter(generateMapButton);
		
		fitShakeMapSourceParam = new BooleanParameter("Fetch Finite Source", true);
		fitShakeMapSourceParam.addParameterChangeListener(this);
		mapParams.addParameter(fitShakeMapSourceParam);
		
		fitMainshockSourceParam = new BooleanParameter("Fit Early Aftershocks", false);
		fitMainshockSourceParam.addParameterChangeListener(this);
		mapParams.addParameter(fitMainshockSourceParam);
		
		publishAdvisoryButton = new ButtonParameter("Publish Advisory", "Publish");
		publishAdvisoryButton.addParameterChangeListener(this);
		publishAdvisoryButton.getEditor().setEditorBorder(BorderFactory.createLineBorder(Color.black, 1));
		publishAdvisoryParams.addParameter(publishAdvisoryButton);
		
		ConsoleWindow console = new ConsoleWindow(true);
		consoleScroll = console.getScrollPane();
		
		JTextArea text = console.getTextArea();
		text.setCaretPosition(0);
		text.setCaretPosition(text.getText().length());

		tabbedPane = new JTabbedPane();
		tabbedPane.setTabLayoutPolicy(JTabbedPane.WRAP_TAB_LAYOUT);
//		tabbedPane.addTab("Console", null, consoleScroll, "View Console"); //moved this to it's own readout panel

		JPanel mainPanel = new JPanel(new BorderLayout());
		JPanel displayPanel = new JPanel(new BorderLayout());
		JPanel paramsPanel = new JPanel(new BorderLayout());
		JPanel dataParamsPanel = new JPanel(new BorderLayout());
		JPanel fitParamsPanel = new JPanel(new BorderLayout());
		JPanel fitParamsSubPanel = new JPanel(new BorderLayout());
		JPanel outputPanel = new JPanel(new BorderLayout());
//		JPanel mapPlotPanel = new JPanel(new BorderLayout());
		
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
		
		
		fitParamsSubPanel.setPreferredSize(new Dimension(paramWidth, (int) (height/10d*(5 + 0.0))));
			mfdEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/10d*(1 + 0.0))));
			fitEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/10d*(4 + 0.0))));
		mapEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/10d*(4 + 0.0))));
		publishAdvisoryEditor.setPreferredSize(new Dimension(paramWidth, (int)(height/10d*(1 + 0.0))));
		
		outputEditor.setPreferredSize(new Dimension((int) (paramWidth/4d), (int) (height/10d * 7)));
		mapPlotEditor.setPreferredSize(new Dimension((int) (paramWidth/4d), (int) (height/10d * 3)));
		
		displayPanel.setPreferredSize(new Dimension(chartWidth, height));
		consoleScroll.setPreferredSize(new Dimension(chartWidth, height));
//		tabbedPane.setPreferredSize(new Dimension(chartWidth, chartHeight));
		paramsPanel.setPreferredSize(new Dimension(2*paramWidth+outputWidth+25, height));
		dataParamsPanel.setPreferredSize(new Dimension(paramWidth, height));
		fitParamsPanel.setPreferredSize(new Dimension(paramWidth, height));
		outputPanel.setPreferredSize(new Dimension(outputWidth+25, height));
//		mapPlotPanel.setPreferredSize(new Dimension(outputWidth+25, height));
		
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
				
		fitParamsSubPanel.add(mfdEditor, BorderLayout.NORTH);
		fitParamsSubPanel.add(fitEditor, BorderLayout.CENTER);
		
//		fitParamsPanel.add(dataEditor, BorderLayout.NORTH);
		fitParamsPanel.add(fitParamsSubPanel, BorderLayout.NORTH);
		fitParamsPanel.add(mapEditor, BorderLayout.CENTER);
		fitParamsPanel.add(publishAdvisoryEditor, BorderLayout.SOUTH);
		
		outputPanel.add(outputEditor, BorderLayout.CENTER);
		outputPanel.add(mapPlotEditor, BorderLayout.SOUTH);

		paramsPanel.add(dataParamsPanel, BorderLayout.WEST);
		paramsPanel.add(fitParamsPanel, BorderLayout.CENTER);
		paramsPanel.add(outputPanel, BorderLayout.EAST);
		
		displayPanel.add(tabbedPane, BorderLayout.CENTER);
		
		//initialize tabs (careful to do it in the same order as in the list)
		tabbedPane.addTab("Console", null, consoleScroll, "Console");
		
		catalogText = new JTextArea();
		catalogText.setEditable(false);
		catalogPane = new JScrollPane(catalogText);
		tabbedPane.addTab("Catalog", null, catalogPane, "Aftershock Catalog");

		tabbedPane.addTab("Epicenters", null, epicenterGraph, "Epicenter Map");
		
		tabbedPane.addTab("Mag/Time Plot", null, magTimeGraph,
				"Aftershock Magnitude vs Time Plot");

		tabbedPane.addTab("Mag/Num Plot", null, magNumGraph,
				"Aftershock Magnitude-Number Distribution");
		
		tabbedPane.addTab("Cml Number Plot", null, cmlNumGraph,
				"Cumulative Number Of Aftershocks with Time");

		tabbedPane.addTab("Model PDFs", null, pdfGraphsPane,
				"Model Probability Distribution Functions");

		tabbedPane.addTab("Forecast MFD", null, forecastMFDPane,
				"Aftershock Magnitude-Number and Probability Distributions");

		tabbedPane.addTab("Forecast Table", null, forecastTablePane,
				"Forecast Table");

		tabbedPane.addTab("Forecast Maps", null, forecastMapPane,
				"Forecast Maps");

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
		
		// resize the window based on how it's fitting on the screen right now.
//		height = mainPanel.getHeight();
//		chartHeight = (int) (0.8*height-30);
//		chartWidth = chartHeight+100;
//		consoleHeight = (int) (0.2*height);

//		setVisible(false);
//		forecastEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/15*5 + 30)));
//		dataEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/15*6 + 30)));
//		mfdEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/15*1 + 30)));
//		fitEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/15*3 + 30)));
//		mapEditor.setPreferredSize(new Dimension(paramWidth, (int) (height/15*2 + 30)));
//		outputEditor.setPreferredSize(new Dimension((int) (paramWidth/4), (int) (height/2)));
//		publishAdvisoryEditor.setPreferredSize(new Dimension(paramWidth, (int)(height/15*1 + 30)));
//		
//		displayPanel.setPreferredSize(new Dimension(chartWidth, height));
//		consoleScroll.setPreferredSize(new Dimension(chartWidth, consoleHeight));
////		tabbedPane.setPreferredSize(new Dimension(chartWidth, chartHeight));
//		paramsPanel.setPreferredSize(new Dimension(2*paramWidth+outputWidth+25, height));
//		dataParamsPanel.setPreferredSize(new Dimension(paramWidth, height));
//		fitParamsPanel.setPreferredSize(new Dimension(paramWidth, height));
//		outputPanel.setPreferredSize(new Dimension(outputWidth+25, height));
		
//		setSize(paramWidth*2 + outputWidth + 25 + chartWidth, height);
//		setVisible(true);
		
		setDefaultCloseOperation(EXIT_ON_CLOSE);
		setTitle("Aftershock Statistics GUI");
		setLocationRelativeTo(null);
		
		workingDir = new File(System.getenv("HOME"));
		accessor = new ComcatAccessor();

		updateRegionParamList(regionTypeParam.getValue(), regionCenterTypeParam.getValue());
		
		setEnableParameterEditing(expertMode);
		refreshTimeWindowEditor();
		setEnableParamsPostFetch(false);
		setEnableParamsPostForecast(false);
			
		printTip(0);
	} //end createAndShowGUI()
	
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
//				regionList.addParameter(regionCenterLocParam);

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
				double radius = 2 + Math.pow((7d/16d)*Math.pow(10, 1.5*mag + 9)/(dS*1e6), 1d/3d) / 300d;
				magSizeFunc.set(i, radius);
			}
			return magSizeFunc;
		}

	// for plotting: how to color the events by magnitude
	private CPT getMagCPT(){
		return getDistCPT();
	}
	
//	private CPT magCPT;
//	private CPT getMagCPT() {
//		if (magCPT != null)
//			return magCPT;
//		EvenlyDiscretizedFunc magSizeFunc = getMagSizeFunc();
//		try {
//			magCPT = GMT_CPT_Files.GMT_WYSIWYG.instance().rescale(
//					magSizeFunc.getMinX(), magSizeFunc.getMaxX());
//		} catch (IOException e) {
//			throw ExceptionUtils.asRuntimeException(e);
//		}
//		return magCPT;
//	}
//	

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

		// These don't display nicely on a terrain map
//		int[][] colorMapValues = new int[][]{
//			{255, 245, 230},
//			{246, 242, 212},
//			{231, 237, 194},
//			{208, 228, 178},
//			{181, 219, 162},
//			{150, 210, 147},
//			{133, 201, 148},
//			{119, 193, 155},
//			{106, 184, 165},
//			{94, 172, 175},
//			{83, 141, 166},
//			{72, 109, 157},
//			{62, 76, 148},
//			{62, 53, 139},
//			{77, 44, 130},
//			{91, 36, 121},
//			{106, 29, 112},
//			{103, 23, 88},
//			{94, 17, 59},
//			{85, 12, 32},
//			{77, 8, 8},
//		};
		
		// try these:
		int[][] colorMapValues = new int[][]{
			{115, 210, 230},
			{109, 196, 227},
			{104, 183, 224},
			{98, 169, 222},
			{92, 155, 219},
			{87, 141, 217},
			{82, 125, 214},
			{77, 109, 212},
			{71, 93, 209},
			{66, 77, 207},
			{61, 61, 204},
			{75, 58, 207},
			{89, 54, 209},
			{103, 50, 212},
			{116, 47, 214},
			{130, 43, 217},
			{150, 39, 212},
			{170, 35, 208},
			{190, 31, 204},
			{210, 27, 199},
			{230, 23, 195},
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
	
//	private CPT getDistCPT() {
//		if (distCPT != null)
//			return distCPT;
//		EvenlyDiscretizedFunc distFunc = getDistFunc();
//		double halfDelta = 0.5*distFunc.getDelta();
//		try {
//			distCPT = new CPT(distFunc.getMinX()-halfDelta, distFunc.getMaxX()+halfDelta, colorMap);
////			distCPT = GMT_CPT_Files.GMT_WYSIWYG.instance().rescale(
////					distFunc.getMinX()-halfDelta, distFunc.getMaxX()+halfDelta);
//		} catch (IOException e) {
//			throw ExceptionUtils.asRuntimeException(e);
//		}
//		return distCPT;
//	}

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

	private static final int tickLabelFontSize = 12;
	private static final int axisLabelFontSize = 14;
	private static final int plotLabelFontSize = 14;

	private static Color[] extra_colors = {Color.GRAY, Color.BLUE, Color.ORANGE, Color.GREEN};

	private static SimpleDateFormat catDateFormat = new SimpleDateFormat("yyyy\tMM\tdd\tHH\tmm\tss");
	private static final TimeZone utc = TimeZone.getTimeZone("UTC");
	static {
		catDateFormat.setTimeZone(utc);
	}

	private static void setupGP(GraphWidget widget) {
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
		//		Preconditions.checkState(eventID != null && !eventID.isEmpty(), "Must supply event ID!");
//		ObsEqkRupture mainshock = null;

		if (eventID == null || eventID.isEmpty())
			System.err.println("Must supply event ID");
		else {
			mainshock = accessor.fetchEvent(eventID);
//			this.mainshock = mainshock;
			if (mainshock == null){
				System.err.println("Event not found: " + eventID);
				setEnableParamsPostFetch(false);
			} else {
//				//fetch details
				ListIterator<?> iter = mainshock.getAddedParametersIterator();
				if(verbose){
					while (iter.hasNext()){
						Parameter<?> param = (Parameter<?>) iter.next();
						System.out.println(param.getName() +":"+ param.getValue());
					}
				}
					
				System.out.println("Mainshock Mag/Lat/Lon/Depth: " + mainshock.getMag() + " " + mainshock.getHypocenterLocation());
//				this.mainshock = mainshock;

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
					aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, region);

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
						aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, region);

						// look again for even larger shocks
						bigAftershocks = aftershocks.getRupsAboveMag(largestShock.getMag() + 0.1);
					}

					// now find centroid 
					if (aftershocks.isEmpty()) {
						if(verbose) System.out.println("No aftershocks found, skipping centroid...");
					} else {
						region = buildRegion(largestShock, getCentroid());
						aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, region);
					}

				} else {
					region = buildRegion(mainshock, null);
					aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, region);
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
	
	private Region buildRegion(ObsEqkRupture event, Location centroid) {
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
			
			return new Region(loc, radius);
		} else  if (type == RegionType.RECTANGULAR) {
			Location lower = new Location(minLatParam.getValue(), minLonParam.getValue());
			Location upper = new Location(maxLatParam.getValue(), maxLonParam.getValue());
			try{
				return new Region(lower, upper);
			} catch (Exception e) {
				System.err.println(e.getMessage());
				return null;
			}
		} else {
			throw new IllegalStateException("Unknown region type: "+type);
		}
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

		// Boolean "changeListenerEnabled" allows the code to suppress paramChangeListener while it builds the model.
//		changeListenerEnabled = false;
		tectonicRegimeParam.setValueAsDefault();
//		changeListenerEnabled = true;
		
		TectonicRegime regime;
		regime = TectonicRegime.GLOBAL_AVERAGE;	//this line is in here just to run faster when debugging (comment out the try/catch if you uncomment this)
		genericFetch = new GenericETAS_ParametersFetch();
		
//		try {
//			if (genericFetch == null)
//				genericFetch = new GenericETAS_ParametersFetch();
//			
//			if(verbose) System.out.println("Determining tectonic regime...");
//			consoleScroll.repaint();
//			regime = genericFetch.getRegion(mainshock.getHypocenterLocation());
//					
//		} catch (RuntimeException e) {
//			System.out.println("Error fetching generic params. Assigning global average values");
//			regime = TectonicRegime.GLOBAL_AVERAGE;	
//						
////			e.printStackTrace();
//		}
		
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
//		cValRangeParam.setValue(new Range(round(genericParams.get_c(),sigDigits), round(genericParams.get_c(),sigDigits))); //c-parameter fixed by default
		
		bParam.setValue(round(genericParams.get_b(),2));
		bParam.getEditor().refreshParamEditor();
		
		magRefParam.setValue(round(genericParams.get_refMag(),2));
		magRefParam.getEditor().refreshParamEditor();
		
		tectonicRegimeParam.setValue(regime);
		tectonicRegimeParam.getEditor().refreshParamEditor();
		
		setMagComplete();
		updateForecastTimes();
		
		// as a courtesy, spit out the decimal days remaining in the origin day
		System.out.println("The mainshock occurred " + String.format("%.4f", getTimeRemainingInUTCDay()) + " days before midnight (UTC)");
		
		if(verbose) System.out.println("Building generic model...");
		// initialize the generic model using an assumed global Mmax of 9.5, max generation depth 100, number of simulations 0
		genericModel = new ETAS_AftershockModel_Generic(mainshock, aftershocks, genericParams, 
				dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
				forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), mcParam.getValue(),
				9.5, 100, 0, fitMSProductivityParam.getValue());

//		changeListenerEnabled = true;
	} //end setupGUI
	
	private void setMagComplete(double mc){
		double b_value = 1; 
		if (ogataMND != null){
			 b_value = ogataMND.calculateBWithFixedMc(mc);
			 if(verbose) System.out.println("mc: " + mc + " b: " + b_value);
		}

		mcParam.setValue(round(mc,2));
		mcParam.getEditor().refreshParamEditor();
		bParam.setValue(round(b_value,2));
		bParam.getEditor().refreshParamEditor();

		link_alpha(); //reset alpha based on b

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
			
			//			mc = ogataMND.get_Mc();
			// iterate to a stable estimate solution
			mc = ogataMND.calculateMcWithFixedB(b_value);
			if(verbose) System.out.println("mc: " + mc + " b: " + b_value);
			
//			int repeat = 0;
//			while(repeat++ < 0){
//				b_value = ogataMND.calculateBWithFixedMc(mc);
//				mc = ogataMND.calculateMcWithFixedB(b_value);
//
//				if(verbose) System.out.println("mc: " + mc + " b: " + b_value);
//			}

		} else {
		 	mc = genericParams.get_refMag() + 0.5;
		}
		
		mcParam.setValue(round(mc,2));
		mcParam.getEditor().refreshParamEditor();
		bParam.setValue(round(b_value,2));
		bParam.getEditor().refreshParamEditor();
		
		link_alpha(); //reset alpha based on b
		
		// this is done in parameter change listener.
//		autoMcParam.setValue(true);
//		computeBButton.setButtonText("Compute b");
//		computeBButton.setInfo("Simultaneously estimate b and Mc.");
//		computeBButton.getEditor().refreshParamEditor();
	}
	
	private void plotAftershockHypocs() {
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		EvenlyDiscretizedFunc magSizeFunc = getMagSizeFunc();
		
		boolean colorByTime = true;
		CPT timeCPT = null;
		
		if (colorByTime) {
			timeCPT = getDistCPT().rescale(0d, dataEndTimeParam.getValue() + 1e-6);
//			try {
//				timeCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, dataEndTimeParam.getValue() + 1e-6); //includ fudge for 0-day forecast
//			} catch (IOException e) {
//				throw ExceptionUtils.asRuntimeException(e);
//			}
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
      		if (subtitle != null)
			epicenterGraph.getGraphPanel().addSubtitle(subtitle);
			
		
	}
	
	private void plotMFDs(IncrementalMagFreqDist mfd, double magComplete) {
		List<XY_DataSet> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		 
//		
//		mfd.setName("Incremental");
//		funcs.add(mfd);
//		chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 6f, Color.BLUE));
//		
		
		
//		//Plot the cumulative density function for Mc (debug)
//		EvenlyDiscretizedFunc mcFunc = ogataMND.getMcCumulativeDensity();
//		if (mcFunc != null){
//			mcFunc.setName("Mc Probability Density");
//			funcs.add(mcFunc);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
//		}
		
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
		
		if(verbose) System.out.println("plot MFD: " + spec.getTitle());
		
		
		if (magNumGraph == null){
			magNumGraph = new GraphWidget(spec);
			Component buttonPanel = magNumGraph.getButtonControlPanel();
			buttonPanel.setVisible(false);
			GraphPanel graphPanel = magNumGraph.getGraphPanel();
			graphPanel.getComponent(2).setVisible(false);
			
			
		} else
			magNumGraph.setPlotSpec(spec);

		magNumGraph.setY_Log(true);

		magNumGraph.setY_AxisRange(plotMinY, plotMaxY);
		magNumGraph.setX_AxisRange(plotMinMag, plotMaxMag);
		setupGP(magNumGraph);


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
		
		if (subtitle != null)
			magTimeGraph.getGraphPanel().addSubtitle(subtitle);
		
		// the color bar disappears if graph is changed with button panel... add some kind of change listener?
	}
	
	private void plotCumulativeNum() {
		double magMin;
		
		if (seqSpecModel != null && timeDepMcParam.getValue() == true)
			magMin = mcParam.getValue();
		else
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
		countFunc.setName("Data: "+(int)countFunc.getMaxY());
		
		double maxY = count + 1;
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(countFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
//		if (seqSpecModel != null && seqSpecModel.nSims != 0 && expertMode ) {
		if (seqSpecModel != null) {
			if (seqSpecModel.get_nSims() != 0 && expertMode) {
				ArbitrarilyDiscretizedFunc expected = getModelCumNumWithLogTimePlot(seqSpecModel, magMin);

				maxY = Math.max(count, expected.getMaxY());

				funcs.add(expected);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, sequence_specific_color));
//				expected.setName("Seq Specific: "+new DecimalFormat("0.#").format(expected.getMaxY()));
				expected.setName("Seq. specific");
				
				ArbitrarilyDiscretizedFunc lower = getFractileCumNumWithLogTimePlot(seqSpecModel, magMin, 0.025);

				funcs.add(lower);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, sequence_specific_color));
				lower.setName("");

				ArbitrarilyDiscretizedFunc upper = getFractileCumNumWithLogTimePlot(seqSpecModel, magMin, 0.975);
				maxY = Math.max(maxY, upper.getMaxY());

				funcs.add(upper);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, sequence_specific_color));
				upper.setName("");

				UncertainArbDiscDataset uncertainFunc = new UncertainArbDiscDataset(expected, lower, upper);
				funcs.add(uncertainFunc);
				uncertainFunc.setName("Seq. specific range");
				PlotLineType plt = PlotLineType.SHADED_UNCERTAIN_TRANS;
//				plt.setAlpha(0.3);
				Color seq_spec_trans_color = new Color(sequence_specific_color.getRed(), sequence_specific_color.getGreen(),sequence_specific_color.getBlue(), (int) (0.3*255) );
				chars.add(new PlotCurveCharacterstics(plt, 1f, seq_spec_trans_color));
			}

			ArbitrarilyDiscretizedFunc fit = getModelFit(seqSpecModel, magMin);
			funcs.add(fit);
			if(expertMode){
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, sequence_specific_color));
				fit.setName("Seq. specific Fit");
			} else {
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, bayesian_color));
				fit.setName("Model Fit");
			}
		} 
		
		if (genericModel != null) {
			if (genericModel.nSims != 0) {
				ArbitrarilyDiscretizedFunc expected = getFractileCumNumWithLogTimePlot(genericModel, magMin, 0.5);
				maxY = Math.max(maxY, expected.getMaxY());

				if (bayesianModel == null || expertMode){
					funcs.add(expected);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, generic_color));
					//				expected.setName("Generic: "+new DecimalFormat("0.#").format(expected.getMaxY()));
					expected.setName("Typical");
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
				PlotLineType plt = PlotLineType.SHADED_UNCERTAIN_TRANS;
//				plt.setAlpha(0.3);
				Color generic_color_trans = new Color(generic_color.getRed(), generic_color.getGreen(),generic_color.getBlue(), (int) (0.3*255) );
				chars.add(new PlotCurveCharacterstics(plt, 1f, generic_color_trans));
			}
			if(expertMode){
				ArbitrarilyDiscretizedFunc fit = getModelFit(genericModel, magMin);
				funcs.add(fit);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, generic_color));
				fit.setName("Typical");
			}
		}
		
		if (bayesianModel != null) {
			if (bayesianModel.nSims != 0) {
				ArbitrarilyDiscretizedFunc expected = getModelCumNumWithLogTimePlot(bayesianModel, magMin);
				maxY = Math.max(maxY, expected.getMaxY());

				funcs.add(expected);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, bayesian_color));
//				expected.setName("Bayesian: "+new DecimalFormat("0.#").format(expected.getMaxY()));
				expected.setName("This forecast");
				
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
				uncertainFunc.setName("Forecast range");
				PlotLineType plt = PlotLineType.SHADED_UNCERTAIN_TRANS;
//				plt.setAlpha(0.3);
				Color bayesian_color_trans = new Color(bayesian_color.getRed(), bayesian_color.getGreen(),bayesian_color.getBlue(), (int) (0.3*255) );
				chars.add(new PlotCurveCharacterstics(plt, 1f, bayesian_color_trans));
			}
			if (expertMode) {
				ArbitrarilyDiscretizedFunc fit = getModelFit(bayesianModel, magMin);
				funcs.add(fit);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, bayesian_color));
				fit.setName("Model fit");
			}
		}

		PlotSpec spec = new PlotSpec(funcs, chars, "Likely number of aftershocks above M"+String.format("%2.1f", magMin), "Days Since Mainshock",
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
	}
	
	private ArbitrarilyDiscretizedFunc getModelFit(ETAS_AftershockModel model, double magMin){
		double tMin = dataStartTimeParam.getValue();
		double tMax = dataEndTimeParam.getValue();
		if (tMax > tMin+1e-9) tMax += -1e-9;
		
		if(verbose) System.out.println("Plotting fit for "+ tMin +" - "+ tMax + " days");
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
			
//			eventIDParam.setName("<custom>");
//			eventIDParam.setValue("Custom");
//			eventIDParam.getEditor().refreshParamEditor();
	
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
			if(verbose) System.out.println("getting generic pdfs");
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
				
				if(verbose) System.out.println("getting bayesian pdfs");
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
		
		if(verbose) System.out.println("getting sequence specific pdfs");
		add1D_PDF(seqSpecModel.getPDF_ams(), "ams-value", amsValExtras);	
		add1D_PDF(seqSpecModel.getPDF_a(), "a-value", aValExtras);
		add1D_PDF(seqSpecModel.getPDF_p(), "p-value", pValExtras);
		add1D_PDF(seqSpecModel.getPDF_logc(), "logc-value", logcValExtras);
		add2D_PDF(seqSpecModel.get2D_PDF_for_a_and_logc(), "a-value", "logc-value");
		add2D_PDF(seqSpecModel.get2D_PDF_for_a_and_p(), "a-value", "p-value");
		add2D_PDF(seqSpecModel.get2D_PDF_for_logc_and_p(), "logc-value", "p-value");
		
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
		
//		Preconditions.checkState(Doubles.isFinite(pdf.getMaxY()), "NaN found in "+pdf.getName());
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
//		try {
//			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(pdf.getMinZ(), pdf.getMaxZ());
//		} catch (IOException e) {
//			throw ExceptionUtils.asRuntimeException(e);
//		}
		
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
//			WindowListener wl = new WindowAdapter() {
//				
//				@Override
//				public void windowClosing(WindowEvent e){
//					stopRequested = true;
//					System.out.println("Calculation interrupted");
//				}
//			};
//			progress.addWindowListener(wl);
			
			pauseRequested = true;
			SwingUtilities.invokeLater(new Runnable() {
			
				@Override
				public void run() {
					if (progress.isVisible())
						progress.setVisible(false);
					
					progress.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
//					progress.setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
					progress.setModalityType(ModalityType.APPLICATION_MODAL);
					
					
					if(steps.length > 1){
//						progress.setIndeterminate(true);
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
//						progress.repaint();
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
//		progress.setIndeterminate(true);
		
		// define common steps
		CalcStep fetchCatalogStep = new CalcStep("Fetching Events",
				"Contacting USGS ComCat Webservice. This is occasionally slow. "
						+ "If it fails, trying again often works.", new Runnable() {

			@Override
			public void run() {
				// set the data window to 0 - now, forecast time to now - one month
//				changeListenerEnabled = false;
				dataStartTimeParam.setValue(0);
				if (dataEndTimeParam.getValue() == null)
					dataEndTimeParam.setValue(dataEndTimeParam.getMax());
				
				fetchEvents();
				if (mainshock == null) stopRequested=true;
//				changeListenerEnabled = true;
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
//				changeListenerEnabled = false;

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
					
					if(seqSpecModel != null)
						seqSpecModel.set_b(b);
					if(bayesianModel != null)
						bayesianModel.set_b(b);
					
					link_alpha(); //reset the alpha parameter, too
					
//					tabbedPane.setSelectedIndex(mag_num_tab_index);
				} else {
					System.err.println("Not enough data to compute b-value");
				}
				
				setEnableParamsPostComputeB(true);
				
//				changeListenerEnabled = true;
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
							genericParams.get_maxMag(), 100, 0, fitMSProductivityParam.getValue());

					resetFitConstraints(genericParams);
				}
			}


		}, true);
		
//		computeAftershockParamsButton
		CalcStep computeAftershockParamStep = new CalcStep("Estimating ETAS parameters", "Computing sequence-specific model. This may take some time...", new Runnable() {

			@Override
			public void run() {
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
								
				if (timeDepMcParam.getValue()) {
					Double rmax = rmaxParam.getValue();
					validateParameter(rmax, "rmax");
				
					//to do: replace with rmax routine from desktop version of the code
					seqSpecModel = new ETAS_AftershockModel_SequenceSpecific(mainshock, aftershocks, rmaxParam.getValue(),
							ETAS_StatsCalc.linspace(amsRange.getLowerBound(), amsRange.getUpperBound(), amsNum), 
							genericParams.get_aSigma(),
							ETAS_StatsCalc.linspace(aRange.getLowerBound(), aRange.getUpperBound(), aNum),
							ETAS_StatsCalc.linspace(pRange.getLowerBound(), pRange.getUpperBound(), pNum),
							ETAS_StatsCalc.logspace(cRange.getLowerBound(), cRange.getUpperBound(), cNum),
							alpha, b, magRefParam.getValue(),	dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
							forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(),	// forecast start time, forecast end time 
							mcParam.getValue(), genericParams.get_maxMag(), 100, 0, fitMSProductivityParam.getValue(), timeDepMcParam.getValue(),
//							new ETAS_AftershockModel_Generic(genericParams), progress);	//max sim mag, max number of generations, number of simulations
							genericModel, progress);	//max sim mag, max number of generations, number of simulations
					
					
				} else {
					seqSpecModel = new ETAS_AftershockModel_SequenceSpecific(mainshock, aftershocks, rmaxParam.getValue(),
							ETAS_StatsCalc.linspace(amsRange.getLowerBound(), amsRange.getUpperBound(), amsNum), 
							genericParams.get_aSigma(),
							ETAS_StatsCalc.linspace(aRange.getLowerBound(), aRange.getUpperBound(), aNum),
							ETAS_StatsCalc.linspace(pRange.getLowerBound(), pRange.getUpperBound(), pNum),
							ETAS_StatsCalc.logspace(cRange.getLowerBound(), cRange.getUpperBound(), cNum),
							alpha, b, magRefParam.getValue(), dataStartTimeParam.getValue(), dataEndTimeParam.getValue(), 
							forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(),
							mcParam.getValue(), genericParams.get_maxMag(), 100, 0, fitMSProductivityParam.getValue(), timeDepMcParam.getValue(),
							null, progress);
							
				}
				System.out.format("Mainshock productivity magnitude: %2.2f\n" , seqSpecModel.getProductivityMag());
			}
			
		}, true);
		
		CalcStep computeBayesStep = new CalcStep("Estimating ETAS parameters", "Computing Bayesian model. This may take some time...", new Runnable() {

			@Override
			public void run() {

				bayesianModel = null;
				if (genericModel != null) {
					bayesianModel = new ETAS_AftershockModel_SequenceSpecific(mainshock, aftershocks, rmaxParam.getValue(),
							seqSpecModel.ams_vec, genericParams.get_aSigma(), seqSpecModel.a_vec, seqSpecModel.p_vec,	seqSpecModel.c_vec,
							seqSpecModel.alpha, seqSpecModel.b, magRefParam.getValue(), dataStartTimeParam.getValue(), dataEndTimeParam.getValue(), 
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
				cValParam.setValue(round(bayesianModel.getMaxLikelihood_c(),6));
				cValParam.getEditor().refreshParamEditor();
			}
		}, false);
		
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
						if(verbose) System.out.println("Computing generic model");
						genericModel.setMagComplete(mcParam.getValue());
						genericModel.computeNewForecast(dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
								forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), 10000);
					}
		}, true);
		
		CalcStep seqSpecForecastStep = new CalcStep("Computing aftershock forecast", "Generating Sequence-Specific forecast. This can take some time...",
				new Runnable() {

					@Override
					public void run() {

						if (seqSpecModel != null){
							seqSpecModel.setMagComplete(mcParam.getValue());
							if(verbose)System.out.println("Sequence Specific model");
							seqSpecModel.computeNewForecast(dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
									forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), 10000);
						}
					}
		}, true);
		
		CalcStep bayesianForecastStep = new CalcStep("Computing aftershock forecast", "Generating Bayesian forecast. This can take some time...",
				new Runnable() {

					@Override
					public void run() {
						
						if (bayesianModel != null){
							bayesianModel.setMagComplete(mcParam.getValue());
							if(verbose)System.out.println("Bayesian model");
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
						System.out.println("Fetching finite source from ShakeMap...");
//						if (!fitMainshockSourceParam.getValue() && faultTrace == null){
						if (fitShakeMapSourceParam.getValue()){
								faultTrace = accessor.fetchShakemapSource(eventIDParam.getValue());
								if(faultTrace == null){
									System.out.println("...ShakeMap finite source not available. Fitting early aftershocks...");
									fitShakeMapSourceParam.setValue(false);
									fitShakeMapSourceParam.getEditor().refreshParamEditor();
									fitMainshockSourceParam.setValue(true);
									fitMainshockSourceParam.getEditor().refreshParamEditor();
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
								if(verbose) System.out.println(faultTrace);
						}
						
						System.out.println("Fetching ShakeMap URL...");
						shakeMapURL = accessor.fetchShakemapURL(eventIDParam.getValue());
						if(D) System.out.println(shakeMapURL);
							
						
					}
		}, true);
		
		CalcStep plotMapStep = new CalcStep("Generating Spatial Forecast", "Generating spatial forecast maps. This can take some time...",
				new Runnable() {

					@Override
					public void run() {
						Stopwatch watch = Stopwatch.createStarted();

						plotRateModel2D(progress);
						
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
//			setChangeListenerEnabled(false); //no cascading updates. All routines must be initiated from EDT, not the changeListener

			if (param == quickForecastButton) {
				System.out.println("Computing quick forecast with defafult settings...");

				CalcStep quickTabRefreshStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
						refreshTabs(forecast_map_tab_index);
					}

				}, true);

				CalcStep quickTabSelectStep = new CalcStep("Plotting","Updating plots...", new Runnable() {

					@Override
					public void run() {
						tabbedPane.setSelectedIndex(forecast_table_tab_index);
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
						computeAftershockParamStep,
						computeBayesStep,
						postParamPlotStep,
						genericForecastStep,
						seqSpecForecastStep,
						bayesianForecastStep,
						postForecastPlotStep,
						plotTableStep,
						fetchSourceStep,
						plotMapStep,
						quickTabRefreshStep,
						quickTabSelectStep,
						tipStep);
				new Thread(run).start();

			} else if (param == eventIDParam || param == dataStartTimeParam || param == dataEndTimeParam
					|| param == regionEditParam) {
				if(verbose) System.out.println("Updating " + param.getName()); 

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
				refreshTimeWindowEditor();

			} else if (param == forecastDurationParam) {
				try{
					if (forecastStartTimeParam.getValue() != null){
						forecastEndTimeParam.setValue(round(forecastStartTimeParam.getValue() + forecastDurationParam.getValue()));
						forecastEndTimeParam.getEditor().refreshParamEditor();
					}
				}finally{

				}
			
			} else if (param == forecastStartTimeParam) {
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
				if(verbose) System.out.println("Updating regionType");
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
				if(verbose) System.out.println("Updating Mc");
				
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

//			} else if (param == autoMcParam) {
//				mcParam.getEditor().setEnabled(!autoMcParam.getValue());
//				mcParam.getEditor().refreshParamEditor();
//			
//				CalcStep setMcStep = new CalcStep("Updating model", "Re-estimating magnitude of completeness...", new Runnable() {
//					@Override
//					public void run() {
//							setMagComplete();
//					}
//				}, true);
//				
//				tabSelectStep = new CalcStep("Plotting","Updating plots...", new Runnable() {
//
//					@Override
//					public void run() {
//						tabbedPane.setSelectedIndex(mag_num_tab_index);
//					}
//
//				}, true);
//				
//				CalcRunnable run;
//				if (autoMcParam.getValue()){
//					run = new CalcRunnable(progress, setMcStep, plotMFDStep, tabSelectStep);
//					new Thread(run).start();
//				} else {
//					boolean prompt;
//				// ask the user if they really want to specify a custom Mc
//					int ret = JOptionPane.showConfirmDialog(this, "Setting the Mc manually may lead to an erroneous forecast,\n" +
//							"and should be done only when the automatic method fails. "
//						+ "\nAre you sure you wish to override the Mc calculation?",
//						"Specify Custom Mc?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
//					prompt = ret == JOptionPane.YES_OPTION;
//					
//					if (!prompt)
//						if (ogataMND != null)
//							mcParam.setValue(ogataMND.get_bValue());
//				}
//					//
////				
////				
//					
			} else if (param == bParam) {
				if(verbose) System.out.println("Updating b-value");

				link_alpha();
//				if (bParam.getValue() != null && bParam.getValue() < 1d)
//					alphaParam.setValue(bParam.getValue()); //link alpha to b if b<1.
//				else
//					alphaParam.setValue(1d);
//
//				if(seqSpecModel != null)
//					seqSpecModel.set_alpha(alphaParam.getValue());
//				if(bayesianModel != null)
//					bayesianModel.set_alpha(alphaParam.getValue());

				CalcStep plotStep = new CalcStep("Plotting", "Building magnitude-frequency distribution plot...", new Runnable() {

					@Override
					public void run() {
						if (tabbedPane.getTabCount() > mag_time_tab_index)
							plotMFDs(aftershockMND, mcParam.getValue());
					}
				}, true);

				CalcRunnable run = new CalcRunnable(progress, plotStep);
				new Thread(run).start();

			} else if (param == fitMSProductivityParam) {
				//grey out the values
				if(fitMSProductivityParam.getValue() != null){
					amsValRangeParam.getEditor().setEnabled(fitMSProductivityParam.getValue());
					amsValNumParam.getEditor().setEnabled(fitMSProductivityParam.getValue());
				}

				setEnableParamsPostAftershockParams(false);

			} else if (param == tectonicRegimeParam) {

				if(verbose) System.out.println("Setting tectonic regime and fetching parameters");
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

				CalcStep updateTectonicRegimeStep = new CalcStep("Updating Model", "Updating generic parameters to " + tectonicRegimeParam.getValue() + "...", new Runnable() {

					@Override
					public void run() {
						genericParams = genericFetch.get(tectonicRegimeParam.getValue());
						Preconditions.checkNotNull(genericParams, "Generic params not found or server error");
						if(verbose) System.out.println("Updating Generic model params to " + tectonicRegimeParam.getValue() + ": "+genericParams);

//						resetFitConstraints(genericParams); //update the constraint editor with the new generic parameters (shift values)
						// update Mc and b
//						double mc = genericParams.get_refMag() + 0.5;
//						mcParam.setValue(mc);
//						mcParam.getEditor().refreshParamEditor();
//						
						setMagComplete();
						
//						genericModel = new ETAS_AftershockModel_Generic(mainshock, aftershocks, genericParams, 
//								dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
//								forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), mcParam.getValue(),
//								9.5, 100, 0, fitMSProductivityParam.getValue());
					}
				}, true);

				CalcRunnable run = new CalcRunnable(progress, updateTectonicRegimeStep, updateModelStep);
				new Thread(run).start();

			} else if (param == amsValRangeParam || param == amsValNumParam) {
				if(verbose) System.out.println("Updating ams-value Range");
				updateRangeParams(amsValRangeParam, amsValNumParam, 51);
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == aValRangeParam || param == aValNumParam) {
				if(verbose) System.out.println("Updating a-value Range");
				updateRangeParams(aValRangeParam, aValNumParam, 51);
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == pValRangeParam || param == pValNumParam) {
				if(verbose) System.out.println("Updating p-value Range");
				updateRangeParams(pValRangeParam, pValNumParam, 51);
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == cValRangeParam || param == cValNumParam) {
				if(verbose) System.out.println("Updating c-value Range");
				updateRangeParams(cValRangeParam, cValNumParam, 51);
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == timeDepMcParam) {
				if(verbose) System.out.println("Time-dependent Mc = " + timeDepMcParam.getValue());
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == rmaxParam){
				if(verbose) System.out.println("Updating rmax-value");
				setEnableParamsPostAftershockParams(false);
				setEnableParamsPostComputeB(false);

			} else if (param == computeAftershockParamsButton) {
				System.out.println("Computing sequence specific fit");

				setEnableParamsPostAftershockParams(false);

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

				CalcRunnable run = new CalcRunnable(progress, computeAftershockParamStep, computeBayesStep, postParamPlotStep, tabRefreshStep, tipStep);
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
						if (seqSpecModel == null)
							printTip(3);
						else
							printTip(5);
					}

				}, true);

				CalcRunnable run = new CalcRunnable(progress, genericForecastStep, seqSpecForecastStep, bayesianForecastStep,
						plotStep, plotTableStep, tabSelectStep, tipStep);
				new Thread(run).start();

			} else if (param == fitMainshockSourceParam) {
				if (fitMainshockSourceParam.getValue()) {
					fitShakeMapSourceParam.setValue(false);
					fitShakeMapSourceParam.getEditor().refreshParamEditor();
				}

			} else if (param == fitShakeMapSourceParam) {
				if (fitShakeMapSourceParam.getValue()){
					fitMainshockSourceParam.setValue(false);
					fitMainshockSourceParam.getEditor().refreshParamEditor();
				}

			} else if (param == generateMapButton) {
				System.out.println("Computing spatial forecast map...");

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

				CalcRunnable run = new CalcRunnable(progress, fetchSourceStep, fetchCitiesStep, plotMapStep, tipStep);
				new Thread(run).start();

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
					plotAllDurationsParam.setValue(true);
					plotAllDurationsParam.getEditor().refreshParamEditor();

					run = new CalcRunnable(progress, postForecastPlotStep, plotMapStep, pubStep, autoSaveStep);	
				} else if (forecastDurationParam.getValue() < 366d) {
					plotAllDurationsParam.setValue(true);
					plotAllDurationsParam.getEditor().refreshParamEditor();
					forecastDurationParam.setValue(366d);
					forecastDurationParam.getEditor().refreshParamEditor();
					updateForecastTimes();
					
					run = new CalcRunnable(progress, 
							computeAftershockParamStep,
							computeBayesStep,
							postParamPlotStep,
							genericForecastStep,
							seqSpecForecastStep,
							bayesianForecastStep,
							postForecastPlotStep,
							plotTableStep,
							fetchSourceStep,
							plotMapStep,
							quickTabRefreshStep,
							pubStep,
							autoSaveStep);
				} else
					run = new CalcRunnable(progress, pubStep, autoSaveStep);
				
				new Thread(run).start();
			} else {
				if(verbose) System.out.println("No action associated with this event.");
			}
			
			setChangeListenerEnabled(true);
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
			
			if (seqSpecModel != null && expertMode) {
				models.add(seqSpecModel);
				names.add("Sequence only");
				colors.add(sequence_specific_color);
			}
			if (genericModel != null) {
				models.add(genericModel);
				names.add("Typical");
				colors.add(generic_color);
			}
			if (bayesianModel != null) {
				models.add(bayesianModel);
				names.add("This forecast");
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

						if(verbose) System.out.println(name + " " + durString);

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

							if(model!=genericModel || expertMode || k != 1){

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

						if(bayesianModel == null || expertMode){
							fractilesFuncs[1].setName(name);
							funcs.add(fractilesFuncs[1]);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
						}

						uncertainFunc.setName(name + " range");
						funcs.add(uncertainFunc);

						PlotLineType plt = PlotLineType.SHADED_UNCERTAIN_TRANS;
						Color c_trans = new Color(c.getRed(), c.getGreen(),c.getBlue(), (int) (0.3*255) );
						PlotCurveCharacterstics uncertainChars = new PlotCurveCharacterstics(plt, 1f, c_trans);
						chars.add(uncertainChars);
					}

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
						if(verbose) System.out.println(yTrack);
					} else {
						yTrack.addValue(0d);
						yTrack.addValue(1d);
						if(verbose) System.out.println(yTrack);
					}

					yTrackProb.addValue(100d);
					if (Double.isFinite(yTrackProb.getMin()) && Double.isFinite(yTrackProb.getMax()) ){
						if(verbose) System.out.println(yTrackProb);
					} else {
						yTrackProb.addValue(0d);
						if(verbose) System.out.println(yTrackProb);
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

					final PlotSpec specNum = new PlotSpec(funcs, chars, "Likely number of aftershocks in the next " + durString, "Magnitude", "Number of aftershocks exceeding magnitude M");
					specNum.setLegendVisible(true);

					final PlotSpec specProb = new PlotSpec(funcsProb, charsProb, "Chance of an aftershock larger than mganitude M in the next " + durString, "Magnitude", "Probability (%)");
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
			if (forecastTablePane == null)
				forecastTablePane = new JTabbedPane();
			else
				while (forecastTablePane.getTabCount() > 0)
					forecastTablePane.removeTabAt(0);
			
			List<ETAS_AftershockModel> models = Lists.newArrayList();
			List<String> names = Lists.newArrayList();
	
			if (bayesianModel != null) {
				models.add(bayesianModel);
				names.add("This forecast");
			}
	
			if (seqSpecModel != null && expertMode){
				models.add(seqSpecModel);
				names.add("Seq. specific");
			}
			
			if (genericModel != null) {
				models.add(genericModel);
				names.add("Typical range");
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
				int maxMag = Math.max(minMag + 3, Math.min(9, (int) Math.ceil(largestShock.getMag()) + 1));
				int nMags = maxMag - minMag + 1;
				if(D) System.out.println("maxMag: " + maxMag + " largestShockMag: " + largestShock.getMag());
				double[] minMags = ETAS_StatsCalc.linspace(minMag, maxMag, nMags);
				
				ETAS_USGS_AftershockForecast forecast = new ETAS_USGS_AftershockForecast(model, minMags, eventDate, startDate, forecastEndTimeParam.getValue());
				JTable jTable = new JTable(forecast.getTableModel());
				jTable.getTableHeader().setFont(jTable.getTableHeader().getFont().deriveFont(Font.BOLD));
				forecastTablePane.addTab(name, jTable);
				
				if(verbose) System.out.println("Took "+watch.elapsed(TimeUnit.SECONDS)+"s to compute aftershock table for "+name);
				watch.stop();
			}
			
			if (tabbedPane.getTabCount() >= forecast_table_tab_index && !tabbedPane.isEnabledAt(forecast_table_tab_index)){
				
				
				tabbedPane.removeTabAt(forecast_table_tab_index);
				tabbedPane.insertTab("Forecast Table", null, forecastTablePane,
						"USGS Forecast Table", forecast_table_tab_index);
			}
			else
				Preconditions.checkState(tabbedPane.getTabCount() > forecast_table_tab_index, "Plots added out of order");
		}

	private void plotRateModel2D(CalcProgressBar progress){
			if (forecastMapPane == null)
				forecastMapPane = new JTabbedPane();
			else
				while (forecastMapPane.getTabCount() > 0)
					forecastMapPane.removeTabAt(0);
	
//			double spacing = 0.0285; 	// grid spacing in degrees
			double spacing = gridSpacingParam.getValue()/111.111; 	// grid spacing in degrees
			double stressDrop = 3.0; 	//MPa
			double mainshockFitDuration = dataEndTimeParam.getValue(); //days
			Double minDays = forecastStartTimeParam.getValue();
			double maxZ = 0;

			// clear out (or initialize) the contourModel
			contourList = new ArrayList<ContourModel>();
			
			// do the full duration models
			String fitType = new String();
			if (fitMainshockSourceParam.getValue())
				fitType = "aftershocks";
			else if (fitShakeMapSourceParam.getValue())
				fitType = "shakemap";
			else
				fitType = "none";
			
			// generate rate model
			GriddedGeoDataSet forecastRateModel;
			if (bayesianModel != null){
				rateModel2D = new ETAS_RateModel2D(bayesianModel);
			} else {
				rateModel2D = new ETAS_RateModel2D(genericModel);
			}
			forecastRateModel = rateModel2D.calculateRateModel(forecastDurationParam.getValue(), spacing, stressDrop, mainshockFitDuration, fitType, faultTrace);
					
			// generate MMI rate model
			GriddedGeoDataSet mmiRateModel;
			int mmiRef = plotMMIParam.getValue();
			System.out.println("MMI level for map: " + mmiRef);
			double poe = plotPOEParam.getValue();
			
			mmiRateModel = rateModel2D.calculateMMIRateModel(mmiRef);
			
			//compute rate sum for checking
			double rateSum = 0;
			for (int j = 0; j < forecastRateModel.size(); j++){
				rateSum += forecastRateModel.get(j);
			}
			if(verbose) System.out.println("rateModel sum: " + rateSum);

			// scoop up the aftershock locations
			PlotSpec oldSpec = epicenterGraph.getPlotSpec();
			List<PlotElement> oldFuncs = (List<PlotElement>) oldSpec.getPlotElems();
			List<PlotCurveCharacterstics> oldChars = oldSpec.getChars();

//			// map city locations for the mapping --TODO where is the labeled XYDataset in OpenSHA?
//			GeoFeatureList cities = loadCities(forecastRateModel.getRegion(), 1000);
//			XY_DataSet cityPlot = new DefaultXY_DataSet();
//			for (GeoFeature city : cities){
//				cityPlot.set(city.loc.getLongitude(), city.loc.getLatitude());
//			}
//			oldFuncs.add(cityPlot);
//			oldChars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 2, Color.BLACK));
			
			// for each duration (day/week/month/year or just what's specified) 
			boolean rateModelPlotted = false;
			for (int i = 0; i < ForecastDuration.values().length; i++){
				// gridded rate model
				GriddedGeoDataSet newForecastRateModel = new GriddedGeoDataSet(forecastRateModel.getRegion(), false);
				for( int j = 0; j < forecastRateModel.size(); j++){
					newForecastRateModel.set(j, forecastRateModel.get(j));
				}
				
				// griddMMImodel
				GriddedGeoDataSet newMMIRateModel = new GriddedGeoDataSet(mmiRateModel.getRegion(), false);
				for( int j = 0; j < mmiRateModel.size(); j++){
					newMMIRateModel.set(j, mmiRateModel.get(j));
				}

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
					//rescale the ratemodel (for speed)
					double targetRate;
					if (bayesianModel == null){
						targetRate = genericModel.getCumNumFractileWithAleatory(new double[]{0.5}, genericModel.magComplete, minDays, maxDays)[0];
					} else {
						targetRate =  bayesianModel.getCumNumFractileWithAleatory(new double[]{0.5}, bayesianModel.magComplete, minDays, maxDays)[0];
					}
					if(verbose) System.out.println("rateModel & targetRate (" + i + "): " + rateSum + " " + targetRate);
					newForecastRateModel.scale(targetRate/rateSum);
					newMMIRateModel.scale(targetRate/rateSum);
					
					// plot the rate only the first time around (should be the longest duration, and all the others are identical except for scaling)
					if (!rateModelPlotted) {
						rateModelPlotted = true;
					
						maxZ = newForecastRateModel.getMaxZ();
						CPT cpt = getProbCPT().rescale(0, maxZ);
						
						List<PolyLine> rateContours = ETAS_RateModel2D.getContours(newForecastRateModel, 20);
						
						List<PlotElement> funcs = Lists.newArrayList();
						List<PlotCurveCharacterstics> chars = Lists.newArrayList();
						for (int j = 1; j < oldFuncs.size() - 1; j++ ){
							if (oldChars.get(j).getSymbol() == PlotSymbol.CIRCLE){
								funcs.add(oldFuncs.get(j));
								chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, (float) (oldChars.get(j).getSymbolWidth()/10), Color.BLACK));
							}
						}
						
						for (int n = 0; n < rateContours.size(); n++){
							XY_DataSet contour = new DefaultXY_DataSet();
							for (PointD pt : rateContours.get(n).PointList)
								contour.set(pt.X,pt.Y);

							funcs.add(contour);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, cpt.getColor((float) rateContours.get(n).Value)));
						}
						Collections.reverse(funcs);
						Collections.reverse(chars);

						PlotSpec spec = new PlotSpec(funcs, chars, "Number of aftershocks larger than M" + String.format("%2.1f", mcParam.getValue())  + " per grid cell in the next " + durString, "Longitude", "Latitude");
						
						GraphWidget forecastMapGraph = new GraphWidget(spec);
						((NumberAxis) forecastMapGraph.getGraphPanel().getXAxis()).setAutoRangeIncludesZero(false);;
						((NumberAxis) forecastMapGraph.getGraphPanel().getYAxis()).setAutoRangeIncludesZero(false);;
						
						Component buttonPanel = forecastMapGraph.getButtonControlPanel();
						buttonPanel.setVisible(false);
						GraphPanel graphPanel = forecastMapGraph.getGraphPanel();
						graphPanel.getComponent(2).setVisible(false);
						
						setupGP(forecastMapGraph);
//						forecastGraph.setBackgroundColor(Color.LIGHT_GRAY);
						
						double regBuff = 0.05;
						Region region = newForecastRateModel.getRegion();
						forecastMapGraph.setAxisRange(region.getMinLon()-regBuff, region.getMaxLon()+regBuff,
								region.getMinLat()-regBuff, region.getMaxLat()+regBuff);
					
						PaintScaleLegend subtitle = XYZGraphPanel.getLegendForCPT(cpt, "Number", axisLabelFontSize, tickLabelFontSize,
					             0, RectangleEdge.RIGHT);
						if (subtitle != null)
							forecastMapGraph.getGraphPanel().addSubtitle(subtitle);

						forecastMapGraph.setName("RateMap (" + durString + ")");
						forecastMapPane.addTab("Rate (" + durString + ")", null, forecastMapGraph);
					}
					
				
//					// draw probGraph
					GriddedGeoDataSet smoothProbModel = rate2PoissProb(newMMIRateModel);
					smoothProbModel.scale(100d);

					CPT cpt2 = getProbCPT().rescale(0, 100);

					double[] contourLevels = ETAS_StatsCalc.linspace(0,100,21);	//specifiy contourLevels directly
					contourLevels[0] = 1;	// no zeros
					contourLevels[contourLevels.length - 1] = 99; 	//no certainty
					
					List<PolyLine> mmiContours = ETAS_RateModel2D.getContours(smoothProbModel, contourLevels);
					
					if(D){
						StringBuilder outString = new StringBuilder();
						for (PolyLine line : mmiContours)
							outString.append(String.format("%2.2f ", line.Value));
						
						System.out.println(outString);
					}
			
//					// write to file? I'd like to put this elsewhere, but need to save all the contours, then...
//					File outFile = new File(workingDir.getAbsolutePath() + "/Contours-"+ durString +".kml");
//					ETAS_RateModel2D.writeContoursAsKML(
//							mmiContours,"Chance of MMI VI (" +durString+")", outFile, cpt2);

					//add contours to the list for later output in publishGraphicalForecast
					contourList.add(new ContourModel(mmiContours,"Chance of exceeding MMI " + decToRoman(mmiRef) + " in the next " +durString, cpt2));
					
					//copy epicenter graph
					List<PlotElement> funcs = Lists.newArrayList();
					List<PlotCurveCharacterstics> chars = Lists.newArrayList();
					for (int j = 1; j < oldFuncs.size() - 1; j++ ){
						if (oldChars.get(j).getSymbol() == PlotSymbol.CIRCLE){
							funcs.add(oldFuncs.get(j));
							chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE,  (float) (oldChars.get(j).getSymbolWidth()/10), Color.BLACK));
						}
					}
										
					for (int n = 0; n < mmiContours.size(); n++){
						XY_DataSet contour = new DefaultXY_DataSet();
						for (PointD pt : mmiContours.get(n).PointList)
							contour.set(pt.X,pt.Y);
						
						funcs.add(contour);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, cpt2.getColor((float) mmiContours.get(n).Value)));
					}
					Collections.reverse(funcs);
					Collections.reverse(chars);

					PlotSpec spec = new PlotSpec(funcs, chars, "Chance of exceeding MMI " + decToRoman(mmiRef) + " in the next " + durString, "Longitude", "Latitude");
					
					GraphWidget forecastGraph = new GraphWidget(spec);
					Component buttonPanel = forecastGraph.getButtonControlPanel();
					buttonPanel.setVisible(false);
					GraphPanel graphPanel = forecastGraph.getGraphPanel();
					graphPanel.getComponent(2).setVisible(false);
					
					double regBuff = 0.05;
					Region region = smoothProbModel.getRegion();
					forecastGraph.setAxisRange(region.getMinLon()-regBuff, region.getMaxLon()+regBuff,
							region.getMinLat()-regBuff, region.getMaxLat()+regBuff);
					
					setupGP(forecastGraph);
					
					PaintScaleLegend subtitle = XYZGraphPanel.getLegendForCPT(cpt2, "Probability (%)", axisLabelFontSize, tickLabelFontSize,
				             0, RectangleEdge.RIGHT);
					if (subtitle != null)
						forecastGraph.getGraphPanel().addSubtitle(subtitle);

					forecastGraph.setName("Probability (" + durString + ")");
					forecastMapPane.addTab("Probability (" + durString + ")", null, forecastGraph);
					
						
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

	private GriddedGeoDataSet rate2PoissProb(GriddedGeoDataSet rateModel){
		double prob;
		for (int i = 0; i < rateModel.size(); i++){
			prob = 1 - Math.exp(-rateModel.get(i));
			rateModel.set(i, prob);
		}
		
		return rateModel;
	}
	
	private String decToRoman(int mmi){
		String str = "";
		switch(mmi) {
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
		workingDirChooser = new JFileChooser(workingDir.getPath() + "/forecast.html");
		
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
			
			if(verbose) System.out.println("Publishing forecast using " + model + "...");
			
			GraphicalForecast graphForecast = new GraphicalForecast(outFile, model, eventDate, startDate);
			graphForecast.setShakeMapURL(shakeMapURL);
			graphForecast.constructForecast();
	
			// print selected figures: Expected number distributions
			for (int i = 0; i < aftershockExpectedNumGraph.size(); i++){
				String file = aftershockExpectedNumGraph.get(i).getName().replaceAll("[^a-zA-Z]",  "").toLowerCase();
				file = outFile.getParent() + "/" + file + ".png"; 
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
	
			// probability of exceeding MMI
			for (int i = 0; i < forecastMapPane.getTabCount(); i++){
				GraphWidget graph = (GraphWidget) forecastMapPane.getComponentAt(i);
				
				String file = graph.getName().replaceAll("[^a-zA-Z]",  "").toLowerCase();
				if (file.contains("probability")){
					file = outFile.getParent() + "/" + file + ".png";
					try {
						System.out.println("Saving forecastProbabilityMap to: " + file);
						graph.saveAsPNG(file);
					} catch (Exception e) {
						System.err.println("Couldn't save forecastProbabilityMap to: " + file);
					}
				} else {
					System.out.println("Skipping " + file);
				}
			}

			// Write Contours to file
			for(ForecastDuration foreDur : ForecastDuration.values()){
				// find the matching contour set
				int nmatch = 0;
				for (ContourModel contours : contourList){
					if (contours.getName().contains(foreDur.toString())){
						nmatch++;
						String name = "contour-" + foreDur.toString();
						File file = new File(outFile.getParent() + "/" + name + ".kml");
						System.out.println("Saving probability contours to: " + file);
						ETAS_RateModel2D.writeContoursAsKML(contours.getContours(), contours.getName(), file, contours.getCPT());
					}
				}
				if (nmatch > 1) System.out.println("More than one set of contours found for duration: " + foreDur);
				else if (nmatch == 0) System.out.println("No contours found for duration: " + foreDur);
			}
			
			
			//write the OFDA logo to the output directory
			// load the data
			String pngFile = "resources/USAID-Logo.png";
			InputStream logoIS = GraphicalForecast.class.getResourceAsStream(pngFile);
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

		if(seqSpecModel != null)
			seqSpecModel.set_alpha(alphaParam.getValue());
		if(bayesianModel != null)
			bayesianModel.set_alpha(alphaParam.getValue());
	}

	private void resetWorkspace() {
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
		seqSpecModel = null;
		faultTrace = null;
		contourList = null;
		
//		if(catalogPane != null)
//			catalogPane.setEnabled(false);
//		if(epicenterGraph != null)
//			epicenterGraph.setEnabled(false);
//		if(magTimeGraph != null)
//			magTimeGraph.setEnabled(false);
//		if(magNumGraph != null)
//			magNumGraph.setEnabled(false);
//		if(cmlNumGraph != null)
//			cmlNumGraph.setEnabled(false);
		
//		if(pdfGraphsPane != null)
//			while (pdfGraphsPane.getTabCount() > 0)
//				pdfGraphsPane.removeTabAt(0);
//		if(forecastMFDPane != null)
//			while (forecastMFDPane.getTabCount() > 0)
//				forecastMFDPane.removeTabAt(0);
//		if(forecastTablePane != null)
//			while (forecastTablePane.getTabCount() > 0)
//				forecastTablePane.removeTabAt(0);
//		if(forecastMapPane != null)
//			while (forecastMapPane.getTabCount() > 0)
//				forecastMapPane.removeTabAt(0);
		
//		for(int i = 0; i < tabbedPane.getTabCount(); i++){
//			tabbedPane.setEnabledAt(i, false);
//		}
		
	}
	
	private void doPostFetchCalculations() {
		// update the grid spacing based on the region size
		gridSpacingParam.setValue(round(Math.sqrt(region.getExtent())/20,1));
		gridSpacingParam.getEditor().refreshParamEditor();
		
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
					"USGS Forecast Table", forecast_table_tab_index);
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
		
//		changeListenerEnabled = false;
		
		dataStartTimeParam.getEditor().setEnabled(!now);
		dataEndTimeParam.getEditor().setEnabled(!now);
		forecastStartTimeParam.getEditor().setEnabled(!now);
		
//		changeListenerEnabled = true;
	}
		
	
	private void resetFitConstraints(GenericETAS_Parameters params){
		if(verbose)System.out.println("Updating fit constraints based on new generic params");
//		changeListenerEnabled = false;
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
		double max_p = pValRangeParam.getValue().getUpperBound();
		
		double min_c = cValRangeParam.getValue().getLowerBound();
		if (min_c <= 1e-6)
			min_c = 1e-6;
		
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
		
		if(restrictParameters && !expertMode && aftershocks.getRupsAboveMag(mc).size() < 10){
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
		
		if(restrictParameters && !expertMode && aftershocks.getRupsAboveMag(mc).size() < 20){
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

		if(restrictParameters && !expertMode && aftershocks.getRupsAboveMag(mc).size() < 30){
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
		if (seqSpecModel != null)
			seqSpecModel = null;
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
		
//		changeListenerEnabled = true;
		
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
		
//		changeListenerEnabled = false;
		
		saveCatalogButton.getEditor().setEnabled(enabled);
		computeBButton.getEditor().setEnabled(enabled);
		tectonicRegimeParam.getEditor().setEnabled(enabled);
		
		// these used to be enabled after computing b but we now allow the user to just use default B
		
		mcParam.getEditor().setEnabled(enabled);
		if (!enabled || (enabled && expertMode)){
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

//		changeListenerEnabled = true;
		
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
		//	 these should not be editable in safeMode, but should be responsive in expertMode	
		amsValParam.getEditor().setEnabled(enabled); // no capability to set directly yet (add in, for custom forecast?)
		aValParam.getEditor().setEnabled(enabled); // no capability to set directly yet
		pValParam.getEditor().setEnabled(enabled); // no capability to set directly yet
		cValParam.getEditor().setEnabled(enabled); // no capability to set directly yet
	}

	private void setEnableParamsPostAftershockParams(boolean enabled) {
//		changeListenerEnabled = false;
		
		if (enabled)
			computeAftershockForecastButton.setButtonText("Run Specific Forecast");
		else
			computeAftershockForecastButton.setButtonText("Run Generic Forecast");
		
		computeAftershockForecastButton.getEditor().refreshParamEditor();
		
		for (int i = aftershock_expected_index; i < tabbedPane.getTabCount() ; i++){
			tabbedPane.setForegroundAt(i, new Color(128,128,128));
			tabbedPane.setEnabledAt(i, false);
		}
		
		setEnableParamsPostForecast(false);
		if (!enabled){
			seqSpecModel = null;
		}
//		changeListenerEnabled = true;
	}
	
	private void setEnableParamsPostForecast(boolean enabled) {
		if (enabled)
			for (int i = forecast_map_tab_index; i < tabbedPane.getTabCount() ; i++){
				tabbedPane.setForegroundAt(i, new Color(128,128,128));
				tabbedPane.setEnabledAt(i, false);
			}
		
		generateMapButton.getEditor().setEnabled(enabled);
		publishAdvisoryButton.getEditor().setEnabled(enabled);
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
		expertMode = false;
			
		//check for arguments
		for (String argument : args){
			System.out.println(argument);
			if (argument.contains("verbose")) verbose = true;
			if (argument.contains("expert")) expertMode = true;
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
		
		System.out.println(formatter.format(origin.getTime()));
		System.out.println(formatter.format(daybreak.getTime()));
		
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

	private GeoFeatureList loadCities(Region region, int minLevel){
		GeoFeatureList cities = new GeoFeatureList();
		System.out.println("Loading Geographic information...");

		// load the data
//		URL citiesURL = GeoFeatureList.class.getResource("resources/worldcities1000.txt"); //world cities with pop >= 1000
//		File cityFile = new File(citiesURL.getFile());

		InputStream citiesIS = GeoFeatureList.class.getResourceAsStream("resources/worldcities1000.txt");
		List<String> lines = new ArrayList<String>();
		try{
//			lines = FileUtils..readLines(cityFile, Charset.defaultCharset());
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
		cities.getFeaturesAboveLevel(minLevel);
		cities.thinFeatures(Math.sqrt(region.getExtent())/10);
		cities.sortByMapLevel(false); //false --> descending order
		if (cities.size() > 10)
			cities.getFeaturesAboveLevel(cities.get(10).mapLevel);
		

		if(verbose && D)
			for(GeoFeature city : cities)
				System.out.println(city);
		
		if(verbose) System.out.println(cities.size() + " cities in mapped region.");
		
		return cities;
				
	}
	
	
	private List<String> tipText;
	private void printTip(int step){
		
		if(tipText == null){
			tipText = new ArrayList<String>();
		}
		
		tipText.add(">> Welcome to the aftershock forecaster. Enter a USGS event ID to get started.");
		tipText.add(">> Specify a forecast start time and duration. Then click \"Fetch Data\" to retrieve the catalog\n  ...or click \"Quick Forecast\" to run the entire forecast automatically with default settings.");
		tipText.add(">> Click \"Compute Model Fit\" to compute sequence-specific model\n  ...or go straight to \"Run Generic Forecast\" to get a generic forecast for this region.");
		tipText.add(">> Click \"Compute Model Fit\" to compute sequence-specific model\n  ...or click \"Render\" to get a generic aftershock rate map.");
		tipText.add(">> Click \"Run Specific Forecast\" to get a sequence-specific forecast.");
		tipText.add(">> Click \"Render\" to generate a sequence-specific aftershock rate map for this forecast.");
		tipText.add(">> Go back and edit the forecast or click \"Publish\" to generate a printable forecast document.");
		
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
