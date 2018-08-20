package scratch.aftershockStatistics;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dialog.ModalityType;
import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Collections;
import java.util.Date;
import java.util.EnumSet;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Map;
import java.util.TimeZone;
import java.util.concurrent.TimeUnit;
import java.util.ArrayDeque;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;

import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.data.Range;
import org.jfree.ui.RectangleEdge;
import org.json.simple.JSONObject;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DatasetBinner;
import scratch.aftershockStatistics.OAFTectonicRegime;
import scratch.aftershockStatistics.USGS_AftershockForecast.Duration;
import scratch.aftershockStatistics.USGS_AftershockForecast.Template;

import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
//import org.opensha.commons.geo.Region;
//import org.opensha.commons.gui.ConsoleWindow;
import org.opensha.commons.gui.plot.GraphWidget;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.param.editor.impl.GriddedParameterListEditor;
import org.opensha.commons.param.editor.impl.ParameterListEditor;
import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.event.ParameterChangeListener;
import org.opensha.commons.param.impl.BooleanParameter;
import org.opensha.commons.param.impl.ButtonParameter;
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.commons.param.impl.EnumParameter;
import org.opensha.commons.param.impl.IntegerParameter;
import org.opensha.commons.param.impl.LocationParameter;
import org.opensha.commons.param.impl.ParameterListParameter;
import org.opensha.commons.param.impl.RangeParameter;
import org.opensha.commons.param.impl.StringParameter;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupListCalc;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
//import org.opensha.sha.gui.infoTools.CalcProgressBar;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import com.lowagie.text.Font;

import gov.usgs.earthquake.product.Product;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
//import scratch.aftershockStatistics.pdl.OAF_Publisher;
import scratch.aftershockStatistics.pdl.PDLProductBuilderOaf;
import scratch.aftershockStatistics.pdl.PDLSender;

import scratch.aftershockStatistics.util.SphLatLon;
import scratch.aftershockStatistics.util.SphRegion;
import scratch.aftershockStatistics.util.GUIConsoleWindow;
import scratch.aftershockStatistics.util.GUICalcStep;
import scratch.aftershockStatistics.util.GUICalcRunnable;
import scratch.aftershockStatistics.util.GUICalcProgressBar;

import scratch.aftershockStatistics.aafs.ServerConfig;
import scratch.aftershockStatistics.aafs.ServerConfigFile;
import scratch.aftershockStatistics.aafs.PDLCmd;

public class AftershockStatsGUI extends JFrame implements ParameterChangeListener {

	// Setting this flag true forces all worker threads to run on the event dispatch thread.

	private boolean forceWorkerEDT = false;

	// Setting this flag true enables the patch for calculation steps that should be on
	// a worker thread, but currently must be on the EDT because they write to the screen.
	// Eventually these calculation steps should be split up.

	private boolean patchWorkerEDT = true;
	
	/*
	 * Data parameters
	 */
	
	private StringParameter eventIDParam;
	private DoubleParameter dataStartTimeParam;
	private DoubleParameter dataEndTimeParam;
	
	private enum RegionType {
		CIRCULAR("Circular"),
		CIRCULAR_WC94("WC 1994 Circular"),
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
		SPECIFIED("Custom Location"),
		CENTROID("Two Step Average Loc");
		
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
	private DoubleParameter radiusParam;
	private DoubleParameter minLatParam;
	private DoubleParameter maxLatParam;
	private DoubleParameter minLonParam;
	private DoubleParameter maxLonParam;
	private DoubleParameter minDepthParam;
	private DoubleParameter maxDepthParam;
	private EnumParameter<RegionCenterType> regionCenterTypeParam;
	private LocationParameter regionCenterLocParam;
	
	private ParameterList regionList;
	private ParameterListParameter regionEditParam;
	
	private ButtonParameter fetchButton;
	
	private JFileChooser loadCatalogChooser;
	private ButtonParameter loadCatalogButton;
	
	private JFileChooser saveCatalogChooser;
	private ButtonParameter saveCatalogButton;
	
	/*
	 * B-value fit parameters
	 */
	
	private DoubleParameter mcParam;
	private DoubleParameter magPrecisionParam;
	private ButtonParameter computeBButton;
	private DoubleParameter bParam;
	
	/*
	 * Aftershock model parameters
	 */
	
	private RangeParameter aValRangeParam;
	private IntegerParameter aValNumParam;
	private RangeParameter pValRangeParam;
	private IntegerParameter pValNumParam;
	private RangeParameter cValRangeParam;
	private IntegerParameter cValNumParam;
	
	private BooleanParameter timeDepMcParam;
	private DoubleParameter gParam;
	private DoubleParameter hParam;
	private DoubleParameter mCatParam;
	
	private ButtonParameter computeAftershockParamsButton;
	
	private DoubleParameter aValParam;
	private DoubleParameter pValParam;
	private DoubleParameter cValParam;
	
	private DoubleParameter forecastStartTimeParam;
	private DoubleParameter forecastEndTimeParam;
	private ButtonParameter forecastStartTimeNowParam;
	
	private ButtonParameter computeAftershockForecastButton;
	
	private JTabbedPane tabbedPane;
	private JScrollPane consoleScroll;
	
	private static final int epicenter_tab_index = 1;
	private static final int mag_num_tab_index = 2;
	private static final int mag_time_tab_index = 3;
	private static final int cml_num_tab_index = 4;
	private static final int catalog_tab_index = 5;
	private static final int pdf_tab_index = 6;
	private static final int aftershock_expected_index = 7;
	private static final int forecast_table_tab_index = 8;
	private GraphWidget epicenterGraph;
	private GraphWidget magNumGraph;
	private GraphWidget magTimeGraph;
	private GraphWidget cmlNumGraph;
	private JTextArea catalogText;
	private JTabbedPane pdfGraphsPane;
	private GraphWidget aftershockExpectedGraph;
	private JTabbedPane forecastTablePane;
	
	private ComcatAccessor accessor;
	private WC1994_MagLengthRelationship wcMagLen;
	
	private SphRegion region;
	private ObsEqkRupture mainshock;
	private ObsEqkRupList aftershocks;
	
	private IncrementalMagFreqDist aftershockMND;
	private double mmaxc;
	
	private RJ_AftershockModel model;
	
	private GenericRJ_ParametersFetch genericFetch = null;
	private GenericRJ_Parameters genericParams = null;
	private RJ_AftershockModel_Generic genericModel = null;
	
	private RJ_AftershockModel_Bayesian bayesianModel = null;
	
	private static final Color generic_color = Color.GREEN.darker();
	private static final Color bayesian_color = Color.RED;
	private static final Color sequence_specific_color = Color.BLUE.darker();
	
	public AftershockStatsGUI() {
		/*
		 * Data parameters
		 */
		ParameterList dataParams = new ParameterList();
		
		eventIDParam = new StringParameter("USGS Event ID");
		eventIDParam.setValue("us20002926");
		eventIDParam.setInfo("Get IDs from https://earthquake.usgs.gov/earthquakes/");
		eventIDParam.addParameterChangeListener(this);
		dataParams.addParameter(eventIDParam);
		
		dataStartTimeParam = new DoubleParameter("Data Start Time", 0d, 36500d, new Double(0d));
		dataStartTimeParam.setUnits("Days");
		dataStartTimeParam.setInfo("Relative to main shock origin time");
		dataStartTimeParam.addParameterChangeListener(this);
		dataParams.addParameter(dataStartTimeParam);
		
		dataEndTimeParam = new DoubleParameter("Data End Time", 0d, 36500d, new Double(7d));
		dataEndTimeParam.setUnits("Days");
		dataEndTimeParam.setInfo("Relative to main shock origin time");
		dataEndTimeParam.addParameterChangeListener(this);
		dataParams.addParameter(dataEndTimeParam);
		
		regionTypeParam = new EnumParameter<AftershockStatsGUI.RegionType>(
				"Region Type", EnumSet.allOf(RegionType.class), RegionType.CIRCULAR_WC94, null);
		regionTypeParam.setInfo("For collecting aftershocks");
		regionTypeParam.addParameterChangeListener(this);
		dataParams.addParameter(regionTypeParam);
		
		// these are inside region editor
		radiusParam = new DoubleParameter("Radius", 0d, 20000d, new Double(20d));
		radiusParam.setUnits("km");
		minLatParam = new DoubleParameter("Min Lat", -90d, 90d, new Double(32d));
		maxLatParam = new DoubleParameter("Max Lat", -90d, 90d, new Double(36d));
		minLonParam = new DoubleParameter("Min Lon", -180d, 180d, new Double(32d));
		maxLonParam = new DoubleParameter("Max Lon", -180d, 180d, new Double(36d));
		minDepthParam = new DoubleParameter("Min Depth", ComcatAccessor.DEFAULT_MIN_DEPTH, ComcatAccessor.DEFAULT_MAX_DEPTH, new Double(ComcatAccessor.DEFAULT_MIN_DEPTH));
		minDepthParam.setUnits("km");
		maxDepthParam = new DoubleParameter("Max Depth", ComcatAccessor.DEFAULT_MIN_DEPTH, ComcatAccessor.DEFAULT_MAX_DEPTH, new Double(ComcatAccessor.DEFAULT_MAX_DEPTH));
		maxDepthParam.setUnits("km");
		regionCenterTypeParam = new EnumParameter<AftershockStatsGUI.RegionCenterType>(
				"Region Center", EnumSet.allOf(RegionCenterType.class), RegionCenterType.CENTROID, null);
		regionCenterLocParam = new LocationParameter("Region Center Location");
		regionCenterTypeParam.addParameterChangeListener(this);
		
		regionList = new ParameterList();
		regionEditParam = new ParameterListParameter("Edit Region", regionList);
		regionEditParam.setInfo("To set more constraints");
		regionEditParam.addParameterChangeListener(this);
		updateRegionParamList(regionTypeParam.getValue(), regionCenterTypeParam.getValue());
		dataParams.addParameter(regionEditParam);
		
		fetchButton = new ButtonParameter("USGS Event Webservice", "Fetch Data");
		fetchButton.setInfo("From USGS ComCat");
		fetchButton.addParameterChangeListener(this);
		dataParams.addParameter(fetchButton);
		
		loadCatalogButton = new ButtonParameter("External Catalog", "Load Catalog");
		loadCatalogButton.setInfo("Load catalog in 10 column format");
		loadCatalogButton.addParameterChangeListener(this);
		dataParams.addParameter(loadCatalogButton);
		
		saveCatalogButton = new ButtonParameter("Aftershock Catalog", "Save Catalog");
		saveCatalogButton.setInfo("Save catalog in 10 column format");
		saveCatalogButton.addParameterChangeListener(this);
		dataParams.addParameter(saveCatalogButton);
		
		mcParam = new DoubleParameter("Mc For Sequence", 0d, 9d);
		mcParam.getConstraint().setNullAllowed(true);
		mcParam.setInfo("Default is Mmaxc+0.5, but user can specify other");
		mcParam.addParameterChangeListener(this);
		dataParams.addParameter(mcParam);
		
		magPrecisionParam = new DoubleParameter("Mag Precision", 0d, 1d, new Double(0.1));
		magPrecisionParam.setInfo("Magnitude rounding applied by network");;
		magPrecisionParam.addParameterChangeListener(this);
		dataParams.addParameter(magPrecisionParam);
		
		computeBButton = new ButtonParameter("Seq. Specific GR b-value", "Compute b (optional)");
		computeBButton.addParameterChangeListener(this);
		dataParams.addParameter(computeBButton);
		
		bParam = new DoubleParameter("b-value", 1d);
		bParam.setValue(null);
		bParam.setInfo("Default is that computed, but user can specify other");
		bParam.addParameterChangeListener(this);
		dataParams.addParameter(bParam);
		
		/*
		 * Fit params
		 */
		
		ParameterList fitParams = new ParameterList();
		
		aValRangeParam = new RangeParameter("a-value range", new Range(-4.5, -0.5));
		aValRangeParam.addParameterChangeListener(this);
		fitParams.addParameter(aValRangeParam);
		
		aValNumParam = new IntegerParameter("a-value num", 1, 10000, new Integer(101));
		aValNumParam.addParameterChangeListener(this);
		fitParams.addParameter(aValNumParam);
		
		pValRangeParam = new RangeParameter("p-value range", new Range(0.9, 2.0));
		pValRangeParam.addParameterChangeListener(this);
		fitParams.addParameter(pValRangeParam);
		
		pValNumParam = new IntegerParameter("p-value num", 1, 10000, new Integer(45));
		pValNumParam.addParameterChangeListener(this);
		fitParams.addParameter(pValNumParam);
		
		cValRangeParam = new RangeParameter("c-value range", new Range(0.018, 0.018));
		cValRangeParam.addParameterChangeListener(this);
		fitParams.addParameter(cValRangeParam);
		
		cValNumParam = new IntegerParameter("c-value num", 1, 10000, new Integer(1));
		cValNumParam.addParameterChangeListener(this);
		fitParams.addParameter(cValNumParam);
		
		timeDepMcParam = new BooleanParameter("Apply time dep. Mc", true);
		timeDepMcParam.addParameterChangeListener(this);
		fitParams.addParameter(timeDepMcParam);
		
		gParam = new DoubleParameter("G", 0.1d, 10d, new Double(0.25));
		gParam.addParameterChangeListener(this);
		fitParams.addParameter(gParam);
		
		hParam = new DoubleParameter("H", 0.25, 2d, new Double(1d));
		hParam.addParameterChangeListener(this);
		fitParams.addParameter(hParam);
		
		mCatParam = new DoubleParameter("Mcat", 1d, 7d, new Double(4.5));
		mCatParam.addParameterChangeListener(this);
		fitParams.addParameter(mCatParam);
		
		computeAftershockParamsButton = new ButtonParameter("Aftershock Params", "Compute");
		computeAftershockParamsButton.addParameterChangeListener(this);
		fitParams.addParameter(computeAftershockParamsButton);
		
		aValParam = new DoubleParameter("a-value", new Double(0d));
		aValParam.setValue(null);
		aValParam.addParameterChangeListener(this);
		fitParams.addParameter(aValParam);
		
		pValParam = new DoubleParameter("p-value", new Double(0d));
		pValParam.setValue(null);
		pValParam.addParameterChangeListener(this);
		fitParams.addParameter(pValParam);
		
		cValParam = new DoubleParameter("c-value", new Double(0d));
		cValParam.setValue(null);
		cValParam.addParameterChangeListener(this);
		fitParams.addParameter(cValParam);
		
		forecastStartTimeNowParam = new ButtonParameter("Set Forecast Start Time", "Set to Now");
		forecastStartTimeNowParam.addParameterChangeListener(this);
		fitParams.addParameter(forecastStartTimeNowParam);
		
		forecastStartTimeParam = new DoubleParameter("Forecast Start Time", 0d, 36500d, new Double(0d));
		forecastStartTimeParam.setUnits("Days");
		forecastStartTimeParam.addParameterChangeListener(this);
		fitParams.addParameter(forecastStartTimeParam);
		
		forecastEndTimeParam = new DoubleParameter("Forecast End Time", 0d, 36500d, new Double(7d));
		forecastEndTimeParam.setUnits("Days");
		forecastEndTimeParam.addParameterChangeListener(this);
		fitParams.addParameter(forecastEndTimeParam);
		
		computeAftershockForecastButton = new ButtonParameter("Aftershock Forecast", "Compute");
		computeAftershockForecastButton.addParameterChangeListener(this);
		fitParams.addParameter(computeAftershockForecastButton);
		
		setEnableParamsPostFetch(false);
		
		ParameterListEditor dataEditor = new ParameterListEditor(dataParams);
		dataEditor.setTitle("Data Parameters");
		
		ParameterListEditor fitEditor = new ParameterListEditor(fitParams);
		fitEditor.setTitle("Aftershock Parameters");
		
		GUIConsoleWindow console = new GUIConsoleWindow(true);
		consoleScroll = console.getScrollPane();
		consoleScroll.setSize(600, 600);
		JTextArea text = console.getTextArea();
		text.setCaretPosition(0);
		text.setCaretPosition(text.getText().length());
		
		tabbedPane = new JTabbedPane();
		tabbedPane.addTab("Console", null, consoleScroll, "View Console");
		
		int paramWidth = 250;
		int chartWidth = 800;
		int height = 900;
		
		JPanel mainPanel = new JPanel(new BorderLayout());
		JPanel paramsPanel = new JPanel(new BorderLayout());
		dataEditor.setPreferredSize(new Dimension(paramWidth, height));
		fitEditor.setPreferredSize(new Dimension(paramWidth, height));
		paramsPanel.add(dataEditor, BorderLayout.WEST);
		paramsPanel.add(fitEditor, BorderLayout.EAST);
		mainPanel.add(paramsPanel, BorderLayout.WEST);
		mainPanel.add(tabbedPane, BorderLayout.CENTER);
		
		setContentPane(mainPanel);
		setSize(250*2+chartWidth, height);
		setDefaultCloseOperation(EXIT_ON_CLOSE);
		setTitle("Aftershock Statistics GUI");
		setLocationRelativeTo(null);
		
		accessor = new ComcatAccessor();
	}
	
	private void updateRegionParamList(RegionType type, RegionCenterType centerType) {
		regionList.clear();
		
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
		
		if (type == RegionType.CIRCULAR || type == RegionType.CIRCULAR_WC94) {
			regionList.addParameter(regionCenterTypeParam);
			if (centerType == RegionCenterType.SPECIFIED)
				regionList.addParameter(regionCenterLocParam);
		}
		
		regionList.addParameter(minDepthParam);
		regionList.addParameter(maxDepthParam);
		
		regionEditParam.getEditor().refreshParamEditor();
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
				radius = wcMagLen.getMedianLength(event.getMag());
				System.out.println("Using Wells & Coppersmith 94 Radius: "+(float)radius+" km");
			}
			
			RegionCenterType centerType = regionCenterTypeParam.getValue();
			Location loc;
			if (centerType == RegionCenterType.EPICENTER)
				loc = event.getHypocenterLocation();
			else if (centerType == RegionCenterType.CENTROID)
				loc = centroid;
			else if (centerType == RegionCenterType.SPECIFIED)
				loc = regionCenterLocParam.getValue();
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
	
	private void fetchEvents() {
		System.out.println("Fetching Events");
		
		String eventID = eventIDParam.getValue();
		Preconditions.checkState(eventID != null && !eventID.isEmpty(), "Must supply event ID!");
		
		mainshock = null;
		ObsEqkRupture mainshock = accessor.fetchEvent(eventID, false, true);	// need extended info for sending to PDL
		Preconditions.checkState(mainshock != null, "Event not found: %s", eventID);
		System.out.println("Mainshock Location: "+mainshock.getHypocenterLocation());
		
		setMainshock(mainshock);
		
		Double minDepth = minDepthParam.getValue();
		validateParameter(minDepth, "min depth");
		Double maxDepth = maxDepthParam.getValue();
		validateParameter(maxDepth, "max depth");
		
		Double minDays = dataStartTimeParam.getValue();
		validateParameter(minDays, "start time");
		Double maxDays = dataEndTimeParam.getValue();
		validateParameter(maxDays, "end time");
		
		// check that end date is before current time
		long eventTime = mainshock.getOriginTime();
		long startTime = eventTime + (long)(minDays*ComcatAccessor.day_millis);
		long endTime = eventTime + (long)(maxDays*ComcatAccessor.day_millis);
		
		Preconditions.checkState(startTime < System.currentTimeMillis(), "Start time is before now!");
		
		if (endTime > System.currentTimeMillis()) {
			double calcMaxDays = (System.currentTimeMillis() - startTime)/ComcatAccessor.day_millis;
			System.out.println("WARNING: End time after current time. Setting max days to: "+calcMaxDays);
			dataEndTimeParam.setValue(calcMaxDays);
			dataEndTimeParam.getEditor().refreshParamEditor();
			maxDays = calcMaxDays;
			validateParameter(maxDays, "end time");
		}
		
		
		if (regionTypeParam.getValue().isCircular()
				&& regionCenterTypeParam.getValue() == RegionCenterType.CENTROID) {
			// first with hypocenter
			region = buildRegion(mainshock, mainshock.getHypocenterLocation());
			
			aftershocks = accessor.fetchAftershocks(mainshock, minDays, maxDays, minDepth, maxDepth, region, region.getPlotWrap());
			
			// now find centroid
			if (aftershocks.isEmpty()) {
				System.out.println("No aftershocks found, skipping centroid");
			} else {
				region = buildRegion(mainshock, getCentroid());
				
				aftershocks = accessor.fetchAftershocks(mainshock, minDays, maxDays, minDepth, maxDepth, region, region.getPlotWrap());
			}
		} else {
			region = buildRegion(mainshock, null);
			
			aftershocks = accessor.fetchAftershocks(mainshock, minDays, maxDays, minDepth, maxDepth, region, region.getPlotWrap());
		}
	}
	
	private void setMainshock(ObsEqkRupture mainshock) {
		this.mainshock = mainshock;
		genericParams = null;
		try (
			ChangeBlock change_block = new ChangeBlock();
		){
			bParam.setValue(null);
		}
		try {
			if (genericFetch == null)
				genericFetch = new GenericRJ_ParametersFetch();
			
			System.out.println("Determining tectonic regime for generic parameters");
			OAFTectonicRegime regime = genericFetch.getRegion(mainshock.getHypocenterLocation());
			Preconditions.checkNotNull(regime, "Regime not found or server error");
			genericParams = genericFetch.get(regime);
			Preconditions.checkNotNull(genericParams, "Generic params not found or server error");
			System.out.println("Generic params for "+regime+": "+genericParams);
			genericModel = new RJ_AftershockModel_Generic(mainshock.getMag(), genericParams);
			// set default values to generic
			pValRangeParam.setValue(new Range(genericParams.get_pValue(), genericParams.get_pValue()));
			pValRangeParam.getEditor().refreshParamEditor();
			cValRangeParam.setValue(new Range(genericParams.get_cValue(), genericParams.get_cValue()));
			cValRangeParam.getEditor().refreshParamEditor();
			try (
				ChangeBlock change_block = new ChangeBlock();
			){
				bParam.setValue(genericModel.get_b());
			}
			bParam.getEditor().refreshParamEditor();
		} catch (RuntimeException e) {
			System.err.println("Error fetching generic params");
			e.printStackTrace();
			genericParams = null;
		}
		// as a courtesy, spit out the decimal days remaining in the origin day
		System.out.println("The mainshock occurred " + String.format("%.4f", getTimeRemainingInUTCDay()) + " days before midnight (UTC)\n");
	}
	
	private Location getCentroid() {
		return AftershockStatsCalc.getSphCentroid(mainshock, aftershocks);
	}
	
	private EvenlyDiscretizedFunc magSizeFunc;
	
	private EvenlyDiscretizedFunc getMagSizeFunc() {
		if (magSizeFunc != null)
			return magSizeFunc;
		
		// size function
		double minMag = 1.25;
		double magDelta = 0.5;
		int numMag = 2*8;
		magSizeFunc = new EvenlyDiscretizedFunc(minMag, numMag, magDelta);
//		double maxMag = magSizeFunc.getMaxX();
//		double minSize = 1d;
//		double maxSize = 20d;
//		double sizeMult = 1.4;
//		double size = minSize;
		
		double dS = 3;
		for (int i=0; i<magSizeFunc.size(); i++) {
			double mag = magSizeFunc.getX(i);
//			double fract = (mag - minMag)/(maxMag - minMag);
//			double size = minSize + fract*(maxSize - minSize);
			
//			magSizeFunc.set(i, size);
//			double radius = Math.pow((7d/16d)*Math.pow(10, 1.5*mag + 9)/(dS*1e6), 1d/3d) / 1000 / 111.111;
			// scale with stress drop, from Nicholas via e-mail 10/26/2015
			double radius = Math.pow((7d/16d)*Math.pow(10, 1.5*mag + 9)/(dS*1e6), 1d/3d) / 300d;
			magSizeFunc.set(i, radius);
//			System.out.println("Mag="+mag+", radius="+radius);
//			size *= sizeMult;
		}
		
		return magSizeFunc;
	}
	
	private EvenlyDiscretizedFunc distFunc;
	
	private EvenlyDiscretizedFunc getDistFunc() {
		if (distFunc == null)
			distFunc = HistogramFunction.getEncompassingHistogram(0d, 199d, 20d);
		
		return distFunc; 
	}
	
	private CPT magCPT;
	
	private CPT getMagCPT() {
		if (magCPT != null)
			return magCPT;
		EvenlyDiscretizedFunc magSizeFunc = getMagSizeFunc();
		try {
			magCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(
					magSizeFunc.getMinX(), magSizeFunc.getMaxX());
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		return magCPT;
	}
	
	private CPT distCPT;
	
	private CPT getDistCPT() {
		if (distCPT != null)
			return distCPT;
		EvenlyDiscretizedFunc distFunc = getDistFunc();
		double halfDelta = 0.5*distFunc.getDelta();
		try {
			distCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(
					distFunc.getMinX()-halfDelta, distFunc.getMaxX()+halfDelta);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		return distCPT;
	}
	
	private static final int tickLabelFontSize = 22;
	private static final int axisLabelFontSize = 24;
	private static final int plotLabelFontSize = 24;
	
	private static void setupGP(GraphWidget widget) {
		widget.setPlotLabelFontSize(plotLabelFontSize);
		widget.setAxisLabelFontSize(axisLabelFontSize);
		widget.setTickLabelFontSize(tickLabelFontSize);
		widget.setBackgroundColor(Color.WHITE);
	}
	
	private void plotAftershockHypocs() {
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		EvenlyDiscretizedFunc magSizeFunc = getMagSizeFunc();
		
		boolean colorByTime = true;
		CPT timeCPT = null;
		
		if (colorByTime) {
			try {
				timeCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, dataEndTimeParam.getValue());
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		PaintScaleLegend subtitle = null;
		
		RuptureSurface mainSurf = mainshock.getRuptureSurface();
		if (mainSurf != null && !mainSurf.isPointSurface()) {
			FaultTrace trace = mainshock.getRuptureSurface().getEvenlyDiscritizedUpperEdge();
			DefaultXY_DataSet traceFunc = new DefaultXY_DataSet();
			traceFunc.setName("Main Shock Trace");
			for(Location loc:trace)
				traceFunc.set(loc.getLongitude(), loc.getLatitude());
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
			EvenlyDiscretizedFunc timeFunc = HistogramFunction.getEncompassingHistogram(0d, timeCPT.getMaxValue()*0.99, 1d);
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
		if (epicenterGraph == null)
			epicenterGraph = new GraphWidget(spec);
		else
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
		
		if (tabbedPane.getTabCount() == epicenter_tab_index)
			tabbedPane.addTab("Epicenters", null, epicenterGraph, "Epicenter Map");
		else
			Preconditions.checkState(tabbedPane.getTabCount() > epicenter_tab_index, "Plots added out of order");
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
	
	private void plotMFDs(IncrementalMagFreqDist mfd, double mmaxc) {
		List<XY_DataSet> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		mfd.setName("Incremental");
		funcs.add(mfd);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 6f, Color.BLUE));
		
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
		cmlMFD.setName("Cumulative");
		funcs.add(cmlMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		boolean addCmlPoints = true;
		
		if (addCmlPoints) {
			ArbitrarilyDiscretizedFunc cmlPoints = new ArbitrarilyDiscretizedFunc();
			cmlPoints.setName(null); // don't show legend
			
			double prevVal = cmlMFD.getY(0);
			
			for (int i=1; i<cmlMFD.size(); i++) {
				double val = cmlMFD.getY(i);
				if (val != prevVal) {
					cmlPoints.set(cmlMFD.getX(i-1), prevVal);
					
					prevVal = val;
				}
			}
			
			funcs.add(cmlPoints);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 6f, Color.BLACK));
		}
		
		double plotMinY = 0.9d;
		double plotMaxY = cmlMFD.getMaxY()+2d;
		
		List<Double> yValsForVerticalLines = Lists.newArrayList(0d, 1e-16, plotMinY, 1d, plotMaxY, 1e3, 2e3, 3e3 ,4e3, 5e3);
		
		// add mainshock mag
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		for (double y : yValsForVerticalLines)
			xy.set(mainshock.getMag(), y);
		xy.setName("Mainshock Mag ("+(float)mainshock.getMag()+")");
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));
		
		// add Mmaxc mag
		xy = new DefaultXY_DataSet();
		for (double y : yValsForVerticalLines)
			xy.set(mmaxc, y);
		xy.setName("Mmaxc ("+(float)mmaxc+")");
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GREEN));
		
		System.out.println("********************Calculating MFD with b: "+bParam.getValue());
		if (bParam.getValue() != null) {
			// add Mc used for b-value calculation
			double mc = mcParam.getValue();
			xy = new DefaultXY_DataSet();
			for (double y : yValsForVerticalLines)
				xy.set(mc, y);
			xy.setName("Mc ("+(float)mc+")");
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.CYAN));
			
			// add best fitting G-R
			double b = bParam.getValue();
			GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(b, 1d, mfd.getMinX(), mfd.getMaxX(), mfd.size());
			// scale to rate at Mc
			int index = mfd.getClosestXIndex(mc);
			gr.scaleToCumRate(index, cmlMFD.getY(index));
//			gr.scaleToIncrRate(index, cmlMFD.getY(index));
			
			gr.setName("G-R b="+(float)b);
//			funcs.add(gr);
			EvenlyDiscretizedFunc cmlGR = gr.getCumRateDistWithOffset();
			funcs.add(cmlGR);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.ORANGE));
			
			plotMaxY = Math.max(plotMaxY, cmlGR.getY(cmlGR.getClosestXIndex(plotMinMag)));
		}
		
		removeEmptyFuncs(funcs, chars);
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Aftershock Mag Num Dist", "Magnitude", "Count");
		spec.setLegendVisible(true);
		
		if (magNumGraph == null)
			magNumGraph = new GraphWidget(spec);
		else
			magNumGraph.setPlotSpec(spec);
		magNumGraph.setY_Log(true);
		
		magNumGraph.setY_AxisRange(plotMinY, plotMaxY);
		magNumGraph.setX_AxisRange(plotMinMag, plotMaxMag);
		setupGP(magNumGraph);
		
		if (tabbedPane.getTabCount() == mag_num_tab_index)
			tabbedPane.addTab("Mag/Num Dist", null, magNumGraph,
					"Aftershock Magnitude vs Number Distribution");
		else
			Preconditions.checkState(tabbedPane.getTabCount() > mag_num_tab_index, "Plots added out of order");
	}
	
	private static void removeEmptyFuncs(List<? extends XY_DataSet> funcs, List<PlotCurveCharacterstics> chars) {
		for (int i=funcs.size(); --i>=0;) {
			if (funcs.get(i).size() == 0) {
				funcs.remove(i);
				chars.remove(i);
			}
		}
	}
	
	private double getTimeSinceMainshock(ObsEqkRupture rup) {
		long ms = mainshock.getOriginTime();
		long as = rup.getOriginTime();
		long delta = as - ms;
		return (double)delta/(1000*60*60*24);
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
			
			for (int i=0; i<aftershocks.size(); i++) {
				ObsEqkRupture aftershock = aftershocks.get(i);
				dists.add(LocationUtils.horzDistanceFast(mainshock.getHypocenterLocation(), aftershock.getHypocenterLocation()));
			}
			
			EvenlyDiscretizedFunc distFunc = getDistFunc();
			
			XY_DataSet[][] binnedFuncs = XY_DatasetBinner.bin2D(points, mags, dists, magSizeFunc, distFunc);
			
			CPT distCPT = getDistCPT();
			
			buildFuncsCharsForBinned2D(binnedFuncs, funcs, chars, distCPT, "dist", distFunc, PlotSymbol.FILLED_CIRCLE);
			
			subtitle = XYZGraphPanel.getLegendForCPT(distCPT, "Distance (km)", axisLabelFontSize, tickLabelFontSize,
					0d, RectangleEdge.RIGHT);
		} else {
			XY_DataSet[] magBinnedFuncs = XY_DatasetBinner.bin(points, mags, magSizeFunc);
			
			buildFuncsCharsForBinned(magBinnedFuncs, funcs, chars, PlotSymbol.CIRCLE);
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Magnitude Vs Time", "Days Since Mainshock", "Magnitude");
		
		if (magTimeGraph == null)
			magTimeGraph = new GraphWidget(spec);
		else
			magTimeGraph.setPlotSpec(spec);
		magTimeGraph.setX_AxisRange(-0.75, dataEndTimeParam.getValue()+0.75);
		magTimeGraph.setY_AxisRange(Math.max(0, magTrack.getMin()-1d), magTrack.getMax()+1d);
		setupGP(magTimeGraph);
		if (subtitle != null)
			magTimeGraph.getGraphPanel().addSubtitle(subtitle);
		
		if (tabbedPane.getTabCount() == mag_time_tab_index)
			tabbedPane.addTab("Mag/Time Plot", null, magTimeGraph,
					"Aftershock Magnitude vs Time Plot");
		else
			Preconditions.checkState(tabbedPane.getTabCount() > mag_time_tab_index, "Plots added out of order");
	}
	
	private void plotCumulativeNum() {
		double magMin;
		
		if (model != null && timeDepMcParam.getValue() == true)
			magMin = mCatParam.getValue();
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
		
		double maxY = count;
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(countFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		if (model != null) {
			EvenlyDiscretizedFunc expected = getModelCumNumWithTimePlot(model, magMin);
			
			maxY = Math.max(count, expected.getMaxY());
			
			funcs.add(expected);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, sequence_specific_color));
			
			expected.setName("Seq Specific: "+new DecimalFormat("0.#").format(expected.getMaxY()));
		}
		
		if (genericModel != null) {
			// calculate generic
			
			EvenlyDiscretizedFunc expected = getModelCumNumWithTimePlot(genericModel, magMin);
			
			maxY = Math.max(count, expected.getMaxY());
			
			funcs.add(expected);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, generic_color));
			
			expected.setName("Generic: "+new DecimalFormat("0.#").format(expected.getMaxY()));
			
			if (bayesianModel != null) {
				expected = getModelCumNumWithTimePlot(bayesianModel, magMin);
				
				maxY = Math.max(count, expected.getMaxY());
				
				funcs.add(expected);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, bayesian_color));
				
				expected.setName("Bayesian: "+new DecimalFormat("0.#").format(expected.getMaxY()));
			}
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Cumulative M≥"+(float)magMin, "Days Since Mainshock",
				"Cumulative Number of Aftershocks");
		spec.setLegendVisible(true);
		
		if (cmlNumGraph == null)
			cmlNumGraph = new GraphWidget(spec);
		else
			cmlNumGraph.setPlotSpec(spec);
		
		setupGP(cmlNumGraph);
//		cmlNumGraph.setX_AxisRange(-0.75, dataEndTimeParam.getValue()+0.75);
//		magTimeGraph.setY_AxisRange(Math.max(0, magTrack.getMin()-1d), magTrack.getMax()+1d);
		cmlNumGraph.setY_AxisRange(0d, maxY*1.1);
		
		if (tabbedPane.getTabCount() == cml_num_tab_index)
			tabbedPane.addTab("Cumulative Num Plot", null, cmlNumGraph,
					"Cumulative Number Of Aftershocks Plot");
		else
			Preconditions.checkState(tabbedPane.getTabCount() > cml_num_tab_index, "Plots added out of order");
	}
	
	private EvenlyDiscretizedFunc getModelCumNumWithTimePlot(RJ_AftershockModel model, double magMin) {
		double tMin = dataStartTimeParam.getValue();
		double tMax = Math.max(dataEndTimeParam.getValue(), forecastEndTimeParam.getValue());
		Preconditions.checkState(tMax > tMin);
		double tDelta = (tMax - tMin)/1000d;
		return model.getModalCumNumEventsWithTime(magMin, tMin, tMax, tDelta);
	}
	
	private static SimpleDateFormat catDateFormat = new SimpleDateFormat("yyyy\tMM\tdd\tHH\tmm\tss");
	private static final TimeZone utc = TimeZone.getTimeZone("UTC");
	static {
		catDateFormat.setTimeZone(utc);
	}
	
	private static String getCatalogLine(ObsEqkRupture rup) {
		StringBuilder sb = new StringBuilder();
		Location hypoLoc = rup.getHypocenterLocation();
		sb.append(catDateFormat.format(rup.getOriginTimeCal().getTime())).append("\t");
		sb.append((float)hypoLoc.getLatitude()).append("\t");
		sb.append((float)hypoLoc.getLongitude()).append("\t");
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
	
	private void plotCatalogText() {
		StringBuilder sb = new StringBuilder();
		sb.append("# Year\tMonth\tDay\tHour\tMinute\tSec\tLat\tLon\tDepth\tMagnitude\n");
		sb.append("# Main Shock:\n");
		sb.append("# ").append(getCatalogLine(mainshock)).append("\n");
		for (ObsEqkRupture rup : aftershocks) {
			sb.append(getCatalogLine(rup)).append("\n");
		}
		if (catalogText == null) {
			Preconditions.checkState(tabbedPane.getTabCount() == catalog_tab_index,  "Plots added out of order");
			catalogText = new JTextArea(sb.toString());
			catalogText.setEditable(false);
			JScrollPane pane = new JScrollPane(catalogText);
			tabbedPane.addTab("Catalog", null, pane, "Aftershock Catalog");
		} else {
			catalogText.setText(sb.toString());
		}
	}
	
	private void loadCatalog(File catalogFile) throws IOException {
		List<String> lines = Files.readLines(catalogFile, Charset.defaultCharset());
		ObsEqkRupList myAftershocks = new ObsEqkRupList();
		ObsEqkRupture myMainshock = null;
		for (int i=0; i<lines.size(); i++) {
			String line = lines.get(i).trim();
			if (line.startsWith("#")) {
				if (line.toLowerCase().startsWith("# main")
						&& i < lines.size()-1 && lines.get(i+1).startsWith("#")) {
					// main shock on next line, starting with a #
					String mainshockLine = lines.get(i+1).substring(1).trim();
					System.out.println("Detected mainshock in file: "+mainshockLine);
					try {
						myMainshock = fromCatalogLine(mainshockLine);
					} catch (Exception e) {
						System.err.println("Error loading mainshock");
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
		
		System.out.println("Loaded "+myAftershocks.size()+" aftershocks from file");
		
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
		if (myMainshock != mainshock) {
			// custom mainshock
			eventIDParam.setName("<custom>");
			eventIDParam.getEditor().refreshParamEditor();
		}
		setMainshock(myMainshock);
		aftershocks = myAftershocks;
		region = null;
		setEnableParamsPostFetch(true);
	}
	
	private void plotPDFs() {
		if (pdfGraphsPane == null)
			pdfGraphsPane = new JTabbedPane();
		else
			while (pdfGraphsPane.getTabCount() > 0)
				pdfGraphsPane.removeTabAt(0);
		
		HistogramFunction[] aValExtras = null;
		if (genericModel != null) {
			HistogramFunction genericA = genericModel.getPDF_a();
			genericA.setName("Generic");
			HistogramFunction bayesianA;
			if (bayesianModel != null) {
				bayesianA = bayesianModel.getPDF_a();
				bayesianA.setName("Bayesian");
				aValExtras = new HistogramFunction[] { genericA, bayesianA };
			} else {
				aValExtras = new HistogramFunction[] { genericA };
			}
		}
		add1D_PDF(model.getPDF_a(), "a-value", aValExtras);
		add1D_PDF(model.getPDF_p(), "p-value");
		add1D_PDF(model.getPDF_c(), "c-value");
		add2D_PDF(model.get2D_PDF_for_a_and_c(), "a-value", "c-value");
		add2D_PDF(model.get2D_PDF_for_a_and_p(), "a-value", "p-value");
		add2D_PDF(model.get2D_PDF_for_c_and_p(), "c-value", "p-value");
		
		if (tabbedPane.getTabCount() == pdf_tab_index)
			tabbedPane.addTab("Model PDFs", null, pdfGraphsPane,
					"Aftershock Model Prob Dist Funcs");
		else
			Preconditions.checkState(tabbedPane.getTabCount() > pdf_tab_index, "Plots added out of order");
	}
	
	private static Color[] extra_colors = {Color.GRAY, Color.BLUE, Color.ORANGE, Color.GREEN};
	
	private void add1D_PDF(HistogramFunction pdf, String name, HistogramFunction... extras) {
		if (pdf == null)
			return;
		
		Preconditions.checkState(Doubles.isFinite(pdf.getMaxY()), "NaN found in "+pdf.getName());
		
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
		pdfGraphsPane.addTab(name, null, widget);
	}
	
	private void add2D_PDF(EvenlyDiscrXYZ_DataSet pdf, String name1, String name2) {
		if (pdf == null)
			return;
		
		String title = "PDF for "+name1+" vs "+name2;
		
		Preconditions.checkState(Doubles.isFinite(pdf.getMaxZ()), "NaN found in "+title);
		
		CPT cpt;
		try {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(pdf.getMinZ(), pdf.getMaxZ());
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		XYZPlotSpec spec = new XYZPlotSpec(pdf, cpt, title, name1, name2, "Density");
		
		XYZGraphPanel xyzGP = new XYZGraphPanel();
		pdfGraphsPane.addTab(name1+" vs "+name2, null, xyzGP);
		double xDelta = pdf.getGridSpacingX();
		double yDelta = pdf.getGridSpacingY();
		xyzGP.drawPlot(spec, false, false,
				new Range(pdf.getMinX()-0.5*xDelta, pdf.getMaxX()+0.5*xDelta),
				new Range(pdf.getMinY()-0.5*yDelta, pdf.getMaxY()+0.5*yDelta));
	}
	
	private void plotExpectedAfershockMFDs(GUICalcProgressBar progress) {
		Double minDays = forecastStartTimeParam.getValue();
		validateParameter(minDays, "start time");
		Double maxDays = forecastEndTimeParam.getValue();
		validateParameter(maxDays, "end time");
		
		double minMag;
		if (mainshock.getMag() < 6)
			minMag = 3d;
		else
			minMag = 3d;
		double maxMag = 9d;
		double deltaMag = 0.1;
		int numMag = (int)((maxMag - minMag)/deltaMag + 1.5);
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		List<RJ_AftershockModel> models = Lists.newArrayList();
		List<String> names = Lists.newArrayList();
		List<Color> colors = Lists.newArrayList();
		
		models.add(model);
		names.add("Seq. Specific");
		colors.add(sequence_specific_color);
		
		if (genericModel != null) {
			models.add(genericModel);
			names.add("Generic");
			colors.add(generic_color);
			
			if (bayesianModel != null) {
				// generate Bayesian model
				models.add(bayesianModel);
				names.add("Bayesian");
				colors.add(bayesian_color);
			}
		}
		
		double[] fractiles = { 0.025, 0.975 };
		
		for (int i=0; i<models.size(); i++) {
			String name = names.get(i);
			if (progress != null)
				progress.updateProgress(i, models.size(), "Calculating "+name+"...");
			RJ_AftershockModel model = models.get(i);
			EvenlyDiscretizedFunc mode = model.getModalCumNumMFD(minMag, maxMag, numMag, minDays, maxDays);
			mode.setName(name+" Mode");
			EvenlyDiscretizedFunc mean = model.getMeanCumNumMFD(minMag, maxMag, numMag, minDays, maxDays);
			mean.setName("Mean");
			Color c = colors.get(i);
			
			funcs.add(mode);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
			
			funcs.add(mean);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, c));
			
			EvenlyDiscretizedFunc[] fractilesFuncs =model.getCumNumMFD_FractileWithAleatoryVariability(
					fractiles, minMag, maxMag, numMag, minDays, maxDays);
			
			for (int j= 0; j<fractiles.length; j++) {
				double f = fractiles[j];
				EvenlyDiscretizedFunc fractile = fractilesFuncs[j];
				fractile.setName("p"+(float)(f*100d)+"%");
				
				funcs.add(fractile);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, c));
			}
		}
		
		if (progress != null) {
			//progress.setIndeterminate(true);
			//progress.setProgressMessage("Plotting...");
			progress.setIndeterminate(true, "Plotting...");
		}
		
		// mainshock mag and Bath's law, use evenly discr functions so that it shows up well at all zoom levels
		double mainshockMag = mainshock.getMag();
		double bathsMag = mainshockMag - 1.2;
		DefaultXY_DataSet mainshockFunc = new DefaultXY_DataSet();
		mainshockFunc.setName("Mainshock M="+(float)mainshockMag);
		DefaultXY_DataSet bathsFunc = new DefaultXY_DataSet();
		bathsFunc.setName("Bath's Law M="+(float)bathsMag);
		
		final MinMaxAveTracker yTrack = new MinMaxAveTracker();
		for (PlotElement elem : funcs) {
			if (elem instanceof XY_DataSet) {
				XY_DataSet xy = (XY_DataSet)elem;
				for (Point2D pt : xy)
					if (pt.getY() > 0)
						yTrack.addValue(pt.getY());
			}
		}
		System.out.println(yTrack);
		EvenlyDiscretizedFunc yVals = new EvenlyDiscretizedFunc(yTrack.getMin(), yTrack.getMax(), 20);
		for (int i=0; i<yVals.size(); i++) {
			double y = yVals.getX(i);
			mainshockFunc.set(mainshockMag, y);
			bathsFunc.set(bathsMag, y);
		}
		funcs.add(mainshockFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
		funcs.add(bathsFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		final PlotSpec spec = new PlotSpec(funcs, chars, "Aftershock Forecast", "Magnitude", "Expected Num \u2265 Mag");
		spec.setLegendVisible(true);
		
		Runnable displayRun = new Runnable() {
			
			@Override
			public void run() {
				if (aftershockExpectedGraph == null)
					aftershockExpectedGraph = new GraphWidget(spec);
				else
					aftershockExpectedGraph.setPlotSpec(spec);
				aftershockExpectedGraph.setY_Log(true);
				aftershockExpectedGraph.setY_AxisRange(new Range(yTrack.getMin(), yTrack.getMax()));
				setupGP(aftershockExpectedGraph);
				
				if (tabbedPane.getTabCount() == aftershock_expected_index)
					tabbedPane.addTab("Forecast", null, aftershockExpectedGraph,
							"Aftershock Expected Frequency Plot");
				else
					Preconditions.checkState(tabbedPane.getTabCount() > aftershock_expected_index, "Plots added out of order");
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
	
	private void plotForecastTable(GUICalcProgressBar progress) {
		if (forecastTablePane == null)
			forecastTablePane = new JTabbedPane();
		else
			while (forecastTablePane.getTabCount() > 0)
				forecastTablePane.removeTabAt(0);
		
		List<RJ_AftershockModel> models = Lists.newArrayList();
		List<String> names = Lists.newArrayList();
		
		models.add(model);
		names.add("Seq. Specific");
		
		if (genericModel != null) {
			models.add(genericModel);
			names.add("Generic");
			
			if (bayesianModel != null) {
				// generate Bayesian model
				models.add(bayesianModel);
				names.add("Bayesian");
			}
		}
		
		GregorianCalendar eventDate = mainshock.getOriginTimeCal();
		GregorianCalendar startDate = new GregorianCalendar();
		Double minDays = forecastStartTimeParam.getValue();
		validateParameter(minDays, "start time");
		double startTime = eventDate.getTime().getTime() + minDays*ProbabilityModelsCalc.MILLISEC_PER_DAY;
		startDate.setTimeInMillis((long)startTime);
		
		for (int i=0; i<models.size(); i++) {
			Stopwatch watch = Stopwatch.createStarted();
			RJ_AftershockModel model = models.get(i);
			String name = names.get(i);
			
			if (progress != null)
				progress.updateProgress(i, models.size(), "Calculating "+name+"...");
			
			USGS_AftershockForecast forecast = new USGS_AftershockForecast(model, aftershocks, eventDate, startDate);
			forecastTablePane.addTab(name, new ForecastTablePanel(forecast));
			System.out.println("Took "+watch.elapsed(TimeUnit.SECONDS)+"s to compute aftershock table for "+name);
			watch.stop();
		}
		if (progress != null)
			progress.updateProgress(models.size(), models.size());
		
		if (tabbedPane.getTabCount() == forecast_table_tab_index)
			tabbedPane.addTab("Forecast Table", null, forecastTablePane,
					"USGS Forecast Table");
		else
			Preconditions.checkState(tabbedPane.getTabCount() > forecast_table_tab_index, "Plots added out of order");
	}
	
	private class ForecastTablePanel extends JPanel implements ParameterChangeListener {
		
		private USGS_AftershockForecast forecast;
		
		private ButtonParameter exportButton;
		private ButtonParameter publishButton;
		private EnumParameter<Duration> advisoryDurationParam;
		private EnumParameter<Template> templateParam;
		private BooleanParameter probAboveMainParam;
		private ButtonParameter injectableTextButton;
		
		private JFileChooser chooser;

		private Product pdl_product;
		private Exception pdl_exception;
		
		public ForecastTablePanel(USGS_AftershockForecast forecast) {
			this.forecast = forecast;
			setLayout(new BorderLayout());
			
			ParameterList params = new ParameterList();
			exportButton = new ButtonParameter("JSON", "Export JSON");
			exportButton.addParameterChangeListener(this);
			params.addParameter(exportButton);
			String publish_forecast = "Publish Forecast (Dry Run)";
			switch ((new ServerConfig()).get_pdl_enable()) {
			case ServerConfigFile.PDLOPT_DEV:
				publish_forecast = "Publish Forecast to PDL-Development";
				break;
			case ServerConfigFile.PDLOPT_PROD:
				publish_forecast = "Publish Forecast to PDL-PRODUCTION";
				break;
			}
			//publishButton = new ButtonParameter("USGS PDL", "Publish Forecast");
			publishButton = new ButtonParameter("USGS PDL", publish_forecast);
			publishButton.addParameterChangeListener(this);
			params.addParameter(publishButton);
			advisoryDurationParam = new EnumParameter<USGS_AftershockForecast.Duration>(
					"Advisory Duration", EnumSet.allOf(Duration.class), forecast.getAdvisoryDuration(), null);
			advisoryDurationParam.addParameterChangeListener(this);
			params.addParameter(advisoryDurationParam);
			templateParam = new EnumParameter<USGS_AftershockForecast.Template>(
					"Template", EnumSet.allOf(Template.class), forecast.getTemplate(), null);
			templateParam.addParameterChangeListener(this);
			params.addParameter(templateParam);
			probAboveMainParam = new BooleanParameter("Include Prob ≥ Main", forecast.isIncludeProbAboveMainshock());
			probAboveMainParam.addParameterChangeListener(this);
			params.addParameter(probAboveMainParam);
			injectableTextButton = new ButtonParameter("Injectable Text", "Set text");
			injectableTextButton.addParameterChangeListener(this);
			params.addParameter(injectableTextButton);
			
			this.add(new GriddedParameterListEditor(params, -1, 2), BorderLayout.NORTH);
			JTable jTable = new JTable(forecast.getTableModel());
			jTable.getTableHeader().setFont(jTable.getTableHeader().getFont().deriveFont(Font.BOLD));
			this.add(jTable, BorderLayout.CENTER);
		}

		@Override
		public void parameterChange(ParameterChangeEvent event) {

			// Ensure we are on the event dispatch thread
			if (!( SwingUtilities.isEventDispatchThread() )) {
				throw new IllegalStateException("AftershockStatsGUI.ForecastTablePanel.parameterChange called while not on the event dispatch thread!");
			}

			if (event.getParameter() == exportButton) {
				if (chooser == null)
					chooser = new JFileChooser();
				int ret = chooser.showSaveDialog(this);
				if (ret == JFileChooser.APPROVE_OPTION) {
					File file = chooser.getSelectedFile();
					JSONObject json = null;
					try {
						json = forecast.buildJSON();
					} catch (Exception e) {
						e.printStackTrace();
						String message = ClassUtils.getClassNameWithoutPackage(e.getClass())+": "+e.getMessage();
						JOptionPane.showMessageDialog(this, message, "Error building JSON", JOptionPane.ERROR_MESSAGE);
					}
					if (json != null) {
						try {
							FileWriter fw = new FileWriter(file);
							fw.write(json.toJSONString());
							fw.close();
						} catch (IOException e) {
							e.printStackTrace();
							String message = ClassUtils.getClassNameWithoutPackage(e.getClass())+": "+e.getMessage();
							JOptionPane.showMessageDialog(this, message, "Error writing JSON", JOptionPane.ERROR_MESSAGE);
						}
					}
				}
			} else if (event.getParameter() == publishButton) {
				String userInput = JOptionPane.showInputDialog(this, "Type \"PDL\" and press OK to publish forecast", "Confirm publication", JOptionPane.PLAIN_MESSAGE);
				if (userInput == null || !(userInput.equals("PDL"))) {
					JOptionPane.showMessageDialog(this, "Canceled: Forecast has NOT been sent to PDL", "Publication canceled", JOptionPane.INFORMATION_MESSAGE);
				} else {
					Product product = null;
					try {
						//product = OAF_Publisher.createProduct(mainshock.getEventId(), forecast);
						JSONObject json = forecast.buildJSON();
						String jsonText = json.toJSONString();
						Map<String, String> eimap = ComcatAccessor.extendedInfoToMap (mainshock, ComcatAccessor.EITMOPT_OMIT_NULL_EMPTY);
						String eventNetwork = eimap.get (ComcatAccessor.PARAM_NAME_NETWORK);
						String eventCode = eimap.get (ComcatAccessor.PARAM_NAME_CODE);
						String eventID = mainshock.getEventId();
						long modifiedTime = 0L;
						boolean isReviewed = true;
						product = PDLProductBuilderOaf.createProduct (eventID, eventNetwork, eventCode, isReviewed, jsonText, modifiedTime);
					} catch (Exception e) {
						e.printStackTrace();
						String message = ClassUtils.getClassNameWithoutPackage(e.getClass())+": "+e.getMessage();
						JOptionPane.showMessageDialog(this, message, "Error building product", JOptionPane.ERROR_MESSAGE);
					}
					if (product != null) {

						//  boolean isSent = false;
						//  try {
						//  	//OAF_Publisher.sendProduct(product);
						//  	PDLSender.signProduct(product);
						//  	PDLSender.sendProduct(product, true);
						//  	isSent = true;
						//  } catch (Exception e) {
						//  	e.printStackTrace();
						//  	String message = ClassUtils.getClassNameWithoutPackage(e.getClass())+": "+e.getMessage();
						//  	JOptionPane.showMessageDialog(this, message, "Error sending product", JOptionPane.ERROR_MESSAGE);
						//  }
						//  if (isSent) {
						//  	JOptionPane.showMessageDialog(this, "Success: Forecast has been successfully sent to PDL", "Publication succeeded", JOptionPane.INFORMATION_MESSAGE);
						//  }

						pdl_product = product;
						pdl_exception = null;
						GUICalcStep pdlSendStep = new GUICalcStep("Sending product to PDL", "...", new Runnable() {
							@Override
							public void run() {
								try {
									//OAF_Publisher.sendProduct(product);
									PDLSender.signProduct(pdl_product);
									PDLSender.sendProduct(pdl_product, true);
								} catch (Exception e) {
									pdl_exception = e;
								}
								SwingUtilities.invokeLater(new Runnable() {
									@Override
									public void run() {
										if (pdl_exception == null) {
										  	JOptionPane.showMessageDialog(ForecastTablePanel.this, "Success: Forecast has been successfully sent to PDL", "Publication succeeded", JOptionPane.INFORMATION_MESSAGE);
										} else {
											pdl_exception.printStackTrace();
											String message = ClassUtils.getClassNameWithoutPackage(pdl_exception.getClass())+": "+pdl_exception.getMessage();
											JOptionPane.showMessageDialog(ForecastTablePanel.this, message, "Error sending product", JOptionPane.ERROR_MESSAGE);
										}
									}
								});
							}
						}, forceWorkerEDT);
						GUICalcRunnable run = new GUICalcRunnable(AftershockStatsGUI.this, pdlSendStep);
						new Thread(run).start();

					}
				}
			} else if (event.getParameter() == advisoryDurationParam) {
				forecast.setAdvisoryDuration(advisoryDurationParam.getValue());
			} else if (event.getParameter() == templateParam) {
				forecast.setTemplate(templateParam.getValue());
			} else if (event.getParameter() == probAboveMainParam) {
				forecast.setIncludeProbAboveMainshock(probAboveMainParam.getValue());
			} else if (event.getParameter() == injectableTextButton) {
				String prevText = forecast.getInjectableText();
				if (prevText == null)
					prevText = "";
				JTextArea area = new JTextArea(prevText);
				Dimension size = new Dimension(300, 200);
				area.setPreferredSize(size);
				area.setMinimumSize(size);
				area.setLineWrap(true);
				area.setWrapStyleWord(true);
				int ret = JOptionPane.showConfirmDialog(this, area, "Set Injectable Text", JOptionPane.OK_CANCEL_OPTION);
				if (ret == JOptionPane.OK_OPTION) {
					String text = area.getText();
					if (text.length() == 0)
						text = null;
					forecast.setInjectableText(text);
				}
			}
		}
	}




	// Mechanism for blocking calls to parameterChange.

	private int change_block_count = 0;

	private class ChangeBlock implements AutoCloseable {

		ChangeBlock () {
			if (!( SwingUtilities.isEventDispatchThread() )) {
				throw new IllegalStateException("AftershockStatsGUI.ChangeBlock called while not on the event dispatch thread!");
			}
			synchronized (AftershockStatsGUI.this) {
				++change_block_count;
			}
		}

		@Override
		public void close() {
			synchronized (AftershockStatsGUI.this) {
				--change_block_count;
			}
			return;
		}
	}




	@Override
	public void parameterChange(ParameterChangeEvent event) {

		// Ensure we are on the event dispatch thread
		if (!( SwingUtilities.isEventDispatchThread() )) {
			throw new IllegalStateException("AftershockStatsGUI.parameterChange called while not on the event dispatch thread!");
		}

		synchronized (this) {
			if (change_block_count != 0) {
				return;
			}
		}

		Parameter<?> param = event.getParameter();
		
		final GUICalcProgressBar progress = new GUICalcProgressBar(this, "", "", false);
		//progress.setIndeterminate(true);
		
		if (param == eventIDParam || param == dataStartTimeParam || param == dataEndTimeParam
				|| param == regionEditParam) {
			setEnableParamsPostFetch(false);
		} else if (param == regionTypeParam || param == regionCenterTypeParam) {
			updateRegionParamList(regionTypeParam.getValue(), regionCenterTypeParam.getValue());
			setEnableParamsPostFetch(false);
		} else if (param == fetchButton) {
			setEnableParamsPostFetch(false);
			GUICalcStep fetchStep = new GUICalcStep("Fetching Events",
					"Contacting USGS ComCat Webservice. This is occasionally slow. "
					+ "If it fails, trying again often works.", new Runnable() {
						
						@Override
						public void run() {
							fetchEvents();
						}
					}, patchWorkerEDT);
			GUICalcStep postFetchPlotStep = new GUICalcStep("Plotting Events/Data", "...", new Runnable() {
						
						@Override
						public void run() {
							doPostFetchPlots();
						}
					}, patchWorkerEDT);
			GUICalcRunnable run = new GUICalcRunnable(progress, fetchStep, postFetchPlotStep);
			new Thread(run).start();
		} else if (param == loadCatalogButton) {
			if (loadCatalogChooser == null)
				loadCatalogChooser = new JFileChooser();
			int ret = loadCatalogChooser.showOpenDialog(this);
			if (ret == JFileChooser.APPROVE_OPTION) {
				final File file = loadCatalogChooser.getSelectedFile();
				GUICalcStep loadStep = new GUICalcStep("Loading events", "...", new Runnable() {
					
					@Override
					public void run() {
						try {
							loadCatalog(file);
						} catch (IOException e) {
							ExceptionUtils.throwAsRuntimeException(e);
						}
					}
				}, patchWorkerEDT);
				GUICalcStep postFetchPlotStep = new GUICalcStep("Plotting Events/Data", "...", new Runnable() {
					
					@Override
					public void run() {
						doPostFetchPlots();
					}
				}, patchWorkerEDT);
				GUICalcRunnable run = new GUICalcRunnable(progress, loadStep, postFetchPlotStep);
				new Thread(run).start();
			}
		} else if (param == saveCatalogButton) {
			if (saveCatalogChooser == null)
				saveCatalogChooser = new JFileChooser();
			int ret = saveCatalogChooser.showSaveDialog(this);
			if (ret == JFileChooser.APPROVE_OPTION) {
				FileWriter fw = null;
				try {
					File file = saveCatalogChooser.getSelectedFile();
					fw = new FileWriter(file);
					fw.write(catalogText.getText());
				} catch (IOException e) {
					e.printStackTrace();
					JOptionPane.showMessageDialog(this, e.getMessage(),
							"Error Saving Catalog", JOptionPane.ERROR_MESSAGE);
				} finally {
					if (fw != null) {
						try {
							fw.close();
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
				}
			}
		} else if (param == mcParam || param == magPrecisionParam) {
			setEnableParamsPostComputeB(false);
			if (genericModel != null)
				bParam.setValue(genericModel.get_b());
			else
				bParam.setValue(null);
			bParam.getEditor().refreshParamEditor();
		} else if (param == computeBButton) {
			setEnableParamsPostComputeB(false);
			GUICalcStep bStep = new GUICalcStep("Computing b", "...", new Runnable() {
				
				@Override
				public void run() {
					Double mc = mcParam.getValue();
					validateParameter(mc, "Mc");
					
					double magPrecision = magPrecisionParam.getValue();
					
					ObsEqkRupList filteredRupList = aftershocks.getRupsAboveMag(mc);
					double b = AftershockStatsCalc.getMaxLikelihood_b_value(filteredRupList, mc, magPrecision);
					System.out.println("Num rups ≥ Mc = "+filteredRupList.size());
					System.out.println("Computed b-value: "+b);
					bParam.setValue(b);
					bParam.getEditor().refreshParamEditor();
					
					setEnableParamsPostComputeB(true);
					tabbedPane.setSelectedIndex(mag_num_tab_index);
				}
			}, patchWorkerEDT);
			GUICalcRunnable run = new GUICalcRunnable(progress, bStep);
			new Thread(run).start();
		} else if (param == bParam) {
			setEnableParamsPostAfershockParams(false);
			if (tabbedPane.getTabCount() > mag_time_tab_index)
				plotMFDs(aftershockMND, mmaxc);
			if (tabbedPane.getTabCount() > cml_num_tab_index)
				plotCumulativeNum();
		} else if (param == aValRangeParam || param == aValNumParam) {
			updateRangeParams(aValRangeParam, aValNumParam, 51);
			setEnableParamsPostAfershockParams(false);
		} else if (param == pValRangeParam || param == pValNumParam) {
			updateRangeParams(pValRangeParam, pValNumParam, 45);
			setEnableParamsPostAfershockParams(false);
		} else if (param == cValRangeParam || param == cValNumParam) {
			updateRangeParams(cValRangeParam, cValNumParam, 45);
			setEnableParamsPostAfershockParams(false);
		} else if (param == timeDepMcParam) {
			setEnableParamsPostAfershockParams(false);
			setEnableParamsPostComputeB(true);
		} else if (param == gParam || param == hParam || param == mCatParam) {
			setEnableParamsPostAfershockParams(false);
		} else if (param == computeAftershockParamsButton) {
			setEnableParamsPostAfershockParams(false);
			GUICalcStep computeStep = new GUICalcStep("Computing Aftershock Params", "...", new Runnable() {

				@Override
				public void run() {
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
					
					if (timeDepMcParam.getValue()) {
						Double g = gParam.getValue();
						validateParameter(g, "G");
						
						Double h = hParam.getValue();
						validateParameter(h, "H");
						
						Double mCat = mCatParam.getValue();
						validateParameter(mCat, "Mcat");
						
						model = new RJ_AftershockModel_SequenceSpecific(mainshock, aftershocks, mCat, g, h, b,
								dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
								aRange.getLowerBound(), aRange.getUpperBound(), aNum,
								pRange.getLowerBound(), pRange.getUpperBound(), pNum,
								cRange.getLowerBound(), cRange.getUpperBound(), cNum);
					} else {
						model = new RJ_AftershockModel_SequenceSpecific(mainshock, aftershocks, mc, b,
								dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
								aRange.getLowerBound(), aRange.getUpperBound(), aNum,
								pRange.getLowerBound(), pRange.getUpperBound(), pNum,
								cRange.getLowerBound(), cRange.getUpperBound(), cNum);
					}
				}
				
			}, patchWorkerEDT);
			GUICalcStep plotStep = new GUICalcStep("Plotting Model PDFs", "...", new Runnable() {
				
				@Override
				public void run() {
					aValParam.setValue(model.getMaxLikelihood_a());
					aValParam.getEditor().refreshParamEditor();
					pValParam.setValue(model.getMaxLikelihood_p());
					pValParam.getEditor().refreshParamEditor();
					cValParam.setValue(model.getMaxLikelihood_c());
					cValParam.getEditor().refreshParamEditor();
					
					bayesianModel = null;
					if (genericModel != null) {
						if (RJ_AftershockModel_Bayesian.areModelsEquivalent(model, genericModel))
							bayesianModel = new RJ_AftershockModel_Bayesian(model, genericModel);
						else
							System.out.println("Could not create Bayesian model as sequence specifc and "
									+ "generic models are not equivalent");
					}
					
					plotPDFs();
					setEnableParamsPostAfershockParams(true);
					plotCumulativeNum();
					tabbedPane.setSelectedIndex(pdf_tab_index);
				}
			}, true);
			GUICalcRunnable run = new GUICalcRunnable(progress, computeStep, plotStep);
			new Thread(run).start();
		} else if (param == forecastStartTimeNowParam) {
			SimpleDateFormat df = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss z");
			GregorianCalendar now = new GregorianCalendar();
			System.out.println("Computing delta from mainshock time ("
					+df.format(mainshock.getOriginTimeCal().getTime())+") to now ("+df.format(now.getTime())+")");
			double delta = USGS_AftershockForecast.getDateDelta(mainshock.getOriginTimeCal(), now);
			System.out.println("Delta: "+delta+" days");
			double prevDiff = forecastEndTimeParam.getValue() - forecastStartTimeParam.getValue();
			if (prevDiff <= 0)
				prevDiff = 7;
			forecastStartTimeParam.setValue(delta);
			forecastEndTimeParam.setValue(delta+prevDiff);
			forecastStartTimeParam.getEditor().refreshParamEditor();
			forecastEndTimeParam.getEditor().refreshParamEditor();
		} else if (param == computeAftershockForecastButton) {
			GUICalcStep plotStep = new GUICalcStep("Computing/Plotting Forecast MFDs", "This can take some time...",
					new Runnable() {

						@Override
						public void run() {
							Stopwatch watch = Stopwatch.createStarted();
							plotExpectedAfershockMFDs(progress);
							watch.stop();
							System.out.println(
									"Took "+watch.elapsed(TimeUnit.SECONDS)+"s to compute/plot aftershock MFDs");
						}
				
			}, patchWorkerEDT);
			GUICalcStep tableStep = new GUICalcStep("Computing Forecast Table", "This can take some time...",
					new Runnable() {

						@Override
						public void run() {
							Stopwatch watch = Stopwatch.createStarted();
							plotForecastTable(progress);
							watch.stop();
							System.out.println(
									"Took "+watch.elapsed(TimeUnit.SECONDS)+"s to compute/plot forecast table");
							tabbedPane.setSelectedIndex(aftershock_expected_index);
						}
				
			}, patchWorkerEDT);
			GUICalcRunnable run = new GUICalcRunnable(progress, plotStep, tableStep);
			new Thread(run).start();
		}
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
		
//		System.out.println(formatter.format(origin.getTime()));
//		System.out.println(formatter.format(daybreak.getTime()));
//		
		daysLeftInDay = 1.0 - ((double)(origin.getTimeInMillis() - daybreak.getTimeInMillis()))/ComcatAccessor.day_millis;
		return daysLeftInDay;
	}
	
	private void doPostFetchPlots() {
		aftershockMND = ObsEqkRupListCalc.getMagNumDist(aftershocks, 1.05, 81, 0.1);
		mmaxc = AftershockStatsCalc.getMmaxC(aftershockMND);
		mcParam.setValue(mmaxc+0.5);
		mcParam.getEditor().refreshParamEditor();
		// plots
		plotAftershockHypocs();
		plotMFDs(aftershockMND, mmaxc);
		plotMagVsTime();
		plotCumulativeNum();
		plotCatalogText();
		
		tabbedPane.setSelectedIndex(epicenter_tab_index);

		setEnableParamsPostFetch(true);
	}
	
	private static void validateParameter(Double value, String name) {
		Preconditions.checkState(value != null, "Must specify "+name);
		Preconditions.checkState(Doubles.isFinite(value), name+" must be finite: %s", value);
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
	
	private void validateRange(Range range, int num, String name) {
		Preconditions.checkState(range != null, "Must supply "+name+" range");
		boolean same = range.getLowerBound() == range.getUpperBound();
		if (same)
			Preconditions.checkState(num == 1, "Num must equal 1 for fixed "+name);
		else
			Preconditions.checkState(num > 1, "Num must be >1 for variable "+name);
	}
	
	/**
	 * disables/enables all parameters that are dependent on the fetch step and beyond
	 */
	private void setEnableParamsPostFetch(boolean enabled) {
		saveCatalogButton.getEditor().setEnabled(enabled);
		mcParam.getEditor().setEnabled(enabled);
		magPrecisionParam.getEditor().setEnabled(enabled);
		computeBButton.getEditor().setEnabled(enabled);
		
		// these used to be enabled after computing b but we now allow the user to just use default B
		bParam.getEditor().setEnabled(enabled);
		aValRangeParam.getEditor().setEnabled(enabled);
		aValNumParam.getEditor().setEnabled(enabled);
		pValRangeParam.getEditor().setEnabled(enabled);
		pValNumParam.getEditor().setEnabled(enabled);
		cValRangeParam.getEditor().setEnabled(enabled);
		cValNumParam.getEditor().setEnabled(enabled);
		computeAftershockParamsButton.getEditor().setEnabled(enabled);
		timeDepMcParam.getEditor().setEnabled(enabled);
		gParam.getEditor().setEnabled(enabled && timeDepMcParam.getValue());
		hParam.getEditor().setEnabled(enabled && timeDepMcParam.getValue());
		mCatParam.getEditor().setEnabled(enabled && timeDepMcParam.getValue());
		
		if (!enabled)
			setEnableParamsPostComputeB(enabled);
	}
	
	/**
	 * disables all parameters that are dependent on the compute b step and beyond
	 */
	private void setEnableParamsPostComputeB(boolean enabled) {
		if (!enabled)
			setEnableParamsPostAfershockParams(enabled);
	}
	
	private void setEnableParamsPostAfershockParams(boolean enabled) {
		aValParam.getEditor().setEnabled(false); // no capability to set in model yet
		pValParam.getEditor().setEnabled(false); // no capability to set in model yet
		cValParam.getEditor().setEnabled(false); // no capability to set in model yet
		forecastStartTimeNowParam.getEditor().setEnabled(enabled);
		forecastStartTimeParam.getEditor().setEnabled(enabled);
		forecastEndTimeParam.getEditor().setEnabled(enabled);
		computeAftershockForecastButton.getEditor().setEnabled(enabled);
		if (!enabled)
			model = null;
	}
	
	public static void main(String[] args) {

		// The GUI accepts command-line arguments for configuring PDL access.
		// Complete details are in PDLCmd.java.
		//
		// To select the PDL destination, use one of these:
		// --pdl=dryrun
		// --pdl=dev
		// --pdl=prod
		// If the --pdl option is not specified, the default is --pdl=dev.
		//
		// To specify a PDL key file, use:
		// --privateKey=PRIVATEKEYFILE
		// If the --privateKey option is not specified, then products are sent unsigned.
		//
		// If you want to delete a product, then in addition to the above you should also include all of these:
		// --delete
		// --code=PRODUCTCODE
		// --eventsource=EVENTNETWORK
		// --eventsourcecode=EVENTCODE
		// The value of --code identifies the product that is to be deleted.  The value of --code is typically an event ID.
		// The values of --eventsource and --eventsourcecode identify the event with which the product is associated;
		// these determine which event page displays the product.
		// When deleting a product, the delete is sent to PDL and the program exits without launching the GUI.
		//
		// If you have a JSON file, then you can send it as a product by including:
		// --update=JSONFILENAME
		// --code=PRODUCTCODE
		// --eventsource=EVENTNETWORK
		// --eventsourcecode=EVENTCODE
		// The value of --code identifies the product that is to be sent.  The value of --code is typically an event ID.
		// The product replaces any prior product that was sent with the same --code.
		// The values of --eventsource and --eventsourcecode identify the event with which the product is associated;
		// these determine which event page displays the product.
		// When sending a product, the product is sent to PDL and the program exits without launching the GUI.

		int lo = 0;
		boolean f_config = true;
		boolean f_send = true;
		int pdl_default = ServerConfigFile.PDLOPT_DEV;

		boolean consumed = PDLCmd.exec_pdl_cmd (args, lo, f_config, f_send, pdl_default);

		if (consumed) {
			return;
		}

		// Run the GUI (Must run on the event dispatch thread!)

		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				new AftershockStatsGUI().setVisible(true);
			}
		});

		return;
	}

}
