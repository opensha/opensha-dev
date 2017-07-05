package scratch.aftershockStatisticsETAS;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Dialog.ModalityType;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Collections;
import java.util.Date;
import java.util.EnumSet;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.TimeZone;
import java.util.concurrent.TimeUnit;

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
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DatasetBinner;
import org.opensha.commons.data.siteData.impl.TectonicRegime;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.ConsoleWindow;
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
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupListCalc;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.CalcProgressBar;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import com.lowagie.text.Font;

import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class AftershockStatsGUI_ETAS extends JFrame implements ParameterChangeListener {
	
	/*
	 * Data parameters
	 */
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
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
	
	private DoubleParameter magRefParam;
	private BooleanParameter useETASParam;
	
	private RangeParameter aValRangeParam;
	private IntegerParameter aValNumParam;
	private RangeParameter pValRangeParam;
	private IntegerParameter pValNumParam;
	private RangeParameter cValRangeParam;
	private IntegerParameter cValNumParam;
	
	private BooleanParameter timeDepMcParam;
	private DoubleParameter rmaxParam;
//	private DoubleParameter gParam;
//	private DoubleParameter hParam;
//	private DoubleParameter mCatParam;

	private ButtonParameter computeAftershockParamsButton;
	
	private DoubleParameter aValParam;
	private DoubleParameter pValParam;
	private DoubleParameter cValParam;
	
	private DoubleParameter forecastStartTimeParam;
	private DoubleParameter forecastEndTimeParam;
	
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
	private static final int forecast_map_tab_index = 9;
	private GraphWidget epicenterGraph;
	private GraphWidget magNumGraph;
	private GraphWidget magTimeGraph;
	private GraphWidget cmlNumGraph;
	private JTextArea catalogText;
	private JTabbedPane pdfGraphsPane;
	private GraphWidget aftershockExpectedGraph;
	private JTabbedPane forecastTablePane;
	private JTabbedPane forecastMapPane;
	
	private ComcatAccessor accessor;
	private WC1994_MagLengthRelationship wcMagLen;
	
	private Region region;
	private ObsEqkRupture mainshock;
	private ObsEqkRupList aftershocks;
	
	private IncrementalMagFreqDist aftershockMND;
	private double mmaxc;
	
//	private RJ_AftershockModel model;
	
//	private GenericRJ_ParametersFetch genericFetch = null;
//	private GenericRJ_Parameters genericParams = null;
//	private RJ_AftershockModel_Generic genericModel = null;
//	private RJ_AftershockModel_Bayesian bayesianModel = null;

//	private ETAS_AftershockModel model;
	
	
	
	private GenericETAS_ParametersFetch genericFetch = null;
	private GenericETAS_Parameters genericParams = null;
	private ETAS_AftershockModel_Generic genericModel = null;
	private ETAS_AftershockModel_SequenceSpecific seqSpecModel = null;
	private ETAS_AftershockModel_Bayesian bayesianModel = null;
	
	private static final Color generic_color = Color.GREEN.darker();
	private static final Color bayesian_color = Color.RED;
	private static final Color sequence_specific_color = Color.BLUE.darker();
	
	public AftershockStatsGUI_ETAS() {
		/*
		 * Data parameters
		 */
		ParameterList dataParams = new ParameterList();
		
		eventIDParam = new StringParameter("USGS Event ID");
		eventIDParam.setValue("us20002926");
		eventIDParam.setInfo("Get IDs from https://earthquake.usgs.gov/earthquakes/");
		eventIDParam.addParameterChangeListener(this);
		dataParams.addParameter(eventIDParam);
		
		dataStartTimeParam = new DoubleParameter("Data Start Time", 0d, 3650, new Double(0d));
		dataStartTimeParam.setUnits("Days");
		dataStartTimeParam.setInfo("Relative to main shock origin time");
		dataStartTimeParam.addParameterChangeListener(this);
		dataParams.addParameter(dataStartTimeParam);
		
		dataEndTimeParam = new DoubleParameter("Data End Time", 0d, 3650, new Double(7d));
		dataEndTimeParam.setUnits("Days");
		dataEndTimeParam.setInfo("Relative to main shock origin time");
		dataEndTimeParam.addParameterChangeListener(this);
		dataParams.addParameter(dataEndTimeParam);
		
		regionTypeParam = new EnumParameter<AftershockStatsGUI_ETAS.RegionType>(
				"Region Type", EnumSet.allOf(RegionType.class), RegionType.CIRCULAR_WC94, null);
		regionTypeParam.setInfo("For collecting aftershocks");
		regionTypeParam.addParameterChangeListener(this);
		dataParams.addParameter(regionTypeParam);
		
		// these are inside region editor
		radiusParam = new DoubleParameter("Radius", 0d, 1000, new Double(20));
		radiusParam.setUnits("km");
		minLatParam = new DoubleParameter("Min Lat", -90d, 90d, new Double(32d));
		maxLatParam = new DoubleParameter("Max Lat", -90d, 90d, new Double(36d));
		minLonParam = new DoubleParameter("Min Lon", -180d, 180d, new Double(32d));
		maxLonParam = new DoubleParameter("Max Lon", -180d, 180d, new Double(36d));
		minDepthParam = new DoubleParameter("Min Depth", 0d, 1000d, new Double(0));
		minDepthParam.setUnits("km");
		maxDepthParam = new DoubleParameter("Max Depth", 0d, 1000d, new Double(1000d));
		maxDepthParam.setUnits("km");
		regionCenterTypeParam = new EnumParameter<AftershockStatsGUI_ETAS.RegionCenterType>(
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
		
		useETASParam = new BooleanParameter("Use ETAS model", true);
		useETASParam.addParameterChangeListener(this);
		fitParams.addParameter(useETASParam);
		
		aValRangeParam = new RangeParameter("a-value range", new Range(-4.0, -1));
		aValRangeParam.addParameterChangeListener(this);
		fitParams.addParameter(aValRangeParam);
		
		aValNumParam = new IntegerParameter("a-value num", 1, 10000, new Integer(101));
		aValNumParam.addParameterChangeListener(this);
		fitParams.addParameter(aValNumParam);
		
		pValRangeParam = new RangeParameter("p-value range", new Range(0.5, 2.0));
		pValRangeParam.addParameterChangeListener(this);
		fitParams.addParameter(pValRangeParam);
		
		pValNumParam = new IntegerParameter("p-value num", 1, 10000, new Integer(31));
		pValNumParam.addParameterChangeListener(this);
		fitParams.addParameter(pValNumParam);
		
		cValRangeParam = new RangeParameter("c-value range", new Range(0.001, 1));
		cValRangeParam.addParameterChangeListener(this);
		fitParams.addParameter(cValRangeParam);
		
		cValNumParam = new IntegerParameter("c-value num", 1, 10000, new Integer(31));
		cValNumParam.addParameterChangeListener(this);
		fitParams.addParameter(cValNumParam);
		
		timeDepMcParam = new BooleanParameter("Apply time dep. Mc", true);
		timeDepMcParam.addParameterChangeListener(this);
		fitParams.addParameter(timeDepMcParam);
		
		rmaxParam = new DoubleParameter("rmax", 1d, 1000d, new Double(200));
		rmaxParam.addParameterChangeListener(this);
		fitParams.addParameter(rmaxParam);
		
//		gParam = new DoubleParameter("G", 0.1d, 10d, new Double(0.25));
//		gParam.addParameterChangeListener(this);
//		fitParams.addParameter(gParam);
//		
//		hParam = new DoubleParameter("H", 0.25, 2d, new Double(1d));
//		hParam.addParameterChangeListener(this);
//		fitParams.addParameter(hParam);
//		
//		mCatParam = new DoubleParameter("Mcat", 1d, 7d, new Double(4.5));
//		mCatParam.addParameterChangeListener(this);
//		fitParams.addParameter(mCatParam);
		
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
		
		magRefParam = new DoubleParameter("Reference Magnitude", new Double(0d));
		
		forecastStartTimeParam = new DoubleParameter("Forecast Start Time", 0d, 3650, new Double(0d));
		forecastStartTimeParam.setUnits("Days");
		forecastStartTimeParam.addParameterChangeListener(this);
		fitParams.addParameter(forecastStartTimeParam);
		
		forecastEndTimeParam = new DoubleParameter("Forecast End Time", 0d, 3650, new Double(7d));
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
		
		ConsoleWindow console = new ConsoleWindow(true);
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
		setVisible(true);
		
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
	
	private Region buildRegion(ObsEqkRupture event, Location centroid) {
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
			
			return new Region(loc, radius);
		} else  if (type == RegionType.RECTANGULAR) {
			Location lower = new Location(minLatParam.getValue(), minLonParam.getValue());
			Location upper = new Location(maxLatParam.getValue(), maxLonParam.getValue());
			return new Region(lower, upper);
		} else {
			throw new IllegalStateException("Unknown region type: "+type);
		}
	}
	
	private void fetchEvents() {
		System.out.println("Fetching Events");
		
		String eventID = eventIDParam.getValue();
		Preconditions.checkState(eventID != null && !eventID.isEmpty(), "Must supply event ID!");
		
		mainshock = null;
		ObsEqkRupture mainshock = accessor.fetchEvent(eventID);
		Preconditions.checkState(mainshock != null, "Event not found: %s", eventID);
		System.out.println("Mainshock Location: "+mainshock.getHypocenterLocation());
		this.mainshock = mainshock;
				
		Double minDepth = minDepthParam.getValue();
		validateParameter(minDepth, "min depth");
		Double maxDepth = maxDepthParam.getValue();
		validateParameter(maxDepth, "max depth");
		
		Double minDays = dataStartTimeParam.getValue();
		validateParameter(minDays, "start time");
		Double maxDays = dataEndTimeParam.getValue();
		
		//check for data end time greater than current time and set to current time if smaller
		double elapsedDays = (double) (System.currentTimeMillis() - mainshock.getOriginTime())/ETAS_StatsCalc.MILLISEC_PER_DAY;
		System.out.println(elapsedDays);
		if( elapsedDays < maxDays){
			maxDays = elapsedDays;
			dataEndTimeParam.setValue(maxDays);
			dataEndTimeParam.getEditor().refreshParamEditor();
		}
		validateParameter(maxDays, "end time");
		
		//update the forecast window based on the data window (default, one year from end of data window)
		forecastStartTimeParam.setValue(Math.ceil(maxDays));
		forecastStartTimeParam.getEditor().refreshParamEditor();
		forecastEndTimeParam.setValue(Math.ceil(maxDays)+365);
		forecastEndTimeParam.getEditor().refreshParamEditor();
		
		if (regionTypeParam.getValue().isCircular()
				&& regionCenterTypeParam.getValue() == RegionCenterType.CENTROID) {
			// first with hypocenter
			region = buildRegion(mainshock, mainshock.getHypocenterLocation());
			
			aftershocks = accessor.fetchAftershocks(mainshock, minDays, maxDays, minDepth, maxDepth, region);
			
			// now find centroid
			if (aftershocks.isEmpty()) {
				System.out.println("No aftershocks found, skipping centroid");
			} else {
				region = buildRegion(mainshock, getCentroid());
				
				aftershocks = accessor.fetchAftershocks(mainshock, minDays, maxDays, minDepth, maxDepth, region);
			}
		} else {
			region = buildRegion(mainshock, null);
			
			aftershocks = accessor.fetchAftershocks(mainshock, minDays, maxDays, minDepth, maxDepth, region);
		}
		
		doPostFetchCalculations();
		
		setMainshock(mainshock);
	}
	
	private void setMainshock(ObsEqkRupture mainshock) {
		this.mainshock = mainshock;
		genericParams = null;
		bParam.setValue(null);
		try {
			if (genericFetch == null)
				genericFetch = new GenericETAS_ParametersFetch();
			
			System.out.println("Determining tectonic regime for generic parameters");
			TectonicRegime regime = genericFetch.getRegion(mainshock.getHypocenterLocation());
//			TectonicRegime regime = TectonicRegime.ANSR_DEEPCON;
			
			if(regime == null){
				genericParams = new GenericETAS_Parameters();
				System.out.println("Generic params for global average " + genericParams);
			}else{
				Preconditions.checkNotNull(regime, "Regime not found or server error");
				genericParams = genericFetch.get(regime);
				Preconditions.checkNotNull(genericParams, "Generic params not found or server error");
				System.out.println("Generic params for "+regime+": "+genericParams);
			}
			
			genericModel = new ETAS_AftershockModel_Generic(mainshock, aftershocks, genericParams, 
					dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
					forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), mcParam.getValue(),
					9.5, 100, 0);
			
			
			
			// set default values to generic
			pValRangeParam.setValue(new Range(genericParams.get_pValue(), genericParams.get_pValue()));
			pValRangeParam.getEditor().refreshParamEditor();
			cValRangeParam.setValue(new Range(genericParams.get_cValue(), genericParams.get_cValue()));
			cValRangeParam.getEditor().refreshParamEditor();
			bParam.setValue(genericModel.get_b());
			bParam.getEditor().refreshParamEditor();
			bParam.setValue(genericModel.get_b());
			bParam.getEditor().refreshParamEditor();
			magRefParam.setValue(genericParams.get_refMag());
			magRefParam.getEditor().refreshParamEditor();
//			mcParam.setValue(genericParams.get_refMag());
//			mcParam.getEditor().refreshParamEditor();
//			
			
			
		} catch (RuntimeException e) {
			System.err.println("Error fetching generic params. Default values assigned");
			e.printStackTrace();
//			genericParams = null;
		}
	}
	
	private Location getCentroid() {
		return ETAS_StatsCalc.getCentroid(mainshock, aftershocks);
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
///		if (seqSpecModel != null && timeDepMcParam.getValue() == true)	//removing Mc(t) method to replace with rmax method
//			magMin = mCatParam.getValue();
//		else
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
		
		if (seqSpecModel != null && seqSpecModel.nSims != 0) {
//			EvenlyDiscretizedFunc expected = getModelCumNumWithTimePlot(seqSpecModel, magMin);
			ArbitrarilyDiscretizedFunc expected = getModelCumNumWithLogTimePlot(seqSpecModel, magMin);
			
			maxY = Math.max(count, expected.getMaxY());
			
			funcs.add(expected);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, sequence_specific_color));
			
			expected.setName("Seq Specific: "+new DecimalFormat("0.#").format(expected.getMaxY()));
		}
		
		if (genericModel != null && genericModel.nSims != 0) {
//			EvenlyDiscretizedFunc expected = getModelCumNumWithTimePlot(genericModel, magMin);
//			ArbitrarilyDiscretizedFunc expected = getModelCumNumWithLogTimePlot(genericModel, magMin);
			ArbitrarilyDiscretizedFunc expected = getFractileCumNumWithLogTimePlot(genericModel, magMin, 0.5);
			
			maxY = Math.max(maxY, expected.getMaxY());
			
			funcs.add(expected);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, generic_color));
			
			expected.setName("Generic: "+new DecimalFormat("0.#").format(expected.getMaxY()));
			
			ArbitrarilyDiscretizedFunc lower = getFractileCumNumWithLogTimePlot(genericModel, magMin, 0.025);
			funcs.add(lower);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, generic_color));

			lower.setName("Generic lower 95% confidence");

			ArbitrarilyDiscretizedFunc upper = getFractileCumNumWithLogTimePlot(genericModel, magMin, 0.975);
			funcs.add(upper);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, generic_color));

			upper.setName("Generic upper 95% confidence");
			
		}
		
		if (bayesianModel != null && bayesianModel.nSims != 0) {
//			EvenlyDiscretizedFunc expected = getModelCumNumWithTimePlot(bayesianModel, magMin);
			ArbitrarilyDiscretizedFunc expected = getModelCumNumWithLogTimePlot(bayesianModel, magMin);

			maxY = Math.max(maxY, expected.getMaxY());

			funcs.add(expected);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, bayesian_color));

			expected.setName("Bayesian: "+new DecimalFormat("0.#").format(expected.getMaxY()));
			
			ArbitrarilyDiscretizedFunc lower = getFractileCumNumWithLogTimePlot(bayesianModel, magMin, 0.025);
			funcs.add(lower);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, bayesian_color));

			lower.setName("Bayesian lower 95% confidence");

			ArbitrarilyDiscretizedFunc upper = getFractileCumNumWithLogTimePlot(bayesianModel, magMin, 0.975);
			funcs.add(upper);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, bayesian_color));

			upper.setName("Bayesian upper 95% confidence");
			

			
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
		cmlNumGraph.setY_AxisRange(1d, maxY*1.1);

		//why is this needed?
		if (tabbedPane.getTabCount() == cml_num_tab_index)
			tabbedPane.addTab("Cumulative Num Plot", null, cmlNumGraph,
					"Cumulative Number Of Aftershocks Plot");
		else
			Preconditions.checkState(tabbedPane.getTabCount() > cml_num_tab_index, "Plots added out of order");
	}
	
//	private EvenlyDiscretizedFunc getModelCumNumWithTimePlot(ETAS_AftershockModel model, double magMin) {
//		System.out.println(model);
//		double tMin = dataStartTimeParam.getValue();
//		double tMax = Math.max(dataEndTimeParam.getValue(), forecastEndTimeParam.getValue());
//		Preconditions.checkState(tMax > tMin);
//		double tDelta = 1;
//		return model.getModalCumNumEventsWithTime(magMin, tMin, tMax, tDelta);
//	}
	
	private ArbitrarilyDiscretizedFunc getModelCumNumWithLogTimePlot(ETAS_AftershockModel model, double magMin) {
		return getFractileCumNumWithLogTimePlot(model, magMin, 0.5);
	}
	
	

	private ArbitrarilyDiscretizedFunc getFractileCumNumWithLogTimePlot(ETAS_AftershockModel model, double magMin, double fractile) {
		System.out.println(model);
		double tMin = forecastStartTimeParam.getValue()+1e-9;
		double tMax = forecastEndTimeParam.getValue();
		Preconditions.checkState(tMax > tMin);
		int numPts = 100;
		return model.getFractileCumNumEventsWithLogTime(magMin, tMin, tMax, numPts, fractile);
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

		aftershocks = myAftershocks;
		
		//update data end time with the maximum time in the catalog
		ObsEqkRupture lastAftershock = aftershocks.get(aftershocks.size()-1);
		System.out.println(lastAftershock.getOriginTime());
		System.out.println(myMainshock.getOriginTime());
		double maxDays = (double) (lastAftershock.getOriginTime() - myMainshock.getOriginTime()) / (double) ETAS_StatsCalc.MILLISEC_PER_DAY;
		dataEndTimeParam.setValue(maxDays);
		dataEndTimeParam.getEditor().refreshParamEditor();
		
		//update the forecast window based on the data window (default, one year from end of data window)
		forecastStartTimeParam.setValue(Math.ceil(maxDays));
		forecastStartTimeParam.getEditor().refreshParamEditor();
		forecastEndTimeParam.setValue(Math.ceil(maxDays)+365);
		forecastEndTimeParam.getEditor().refreshParamEditor();
	

		region = null;
		doPostFetchCalculations();
		
		setMainshock(myMainshock);
		
		setEnableParamsPostFetch(true);
		
	}
	
	private void plotPDFs() {
		if (pdfGraphsPane == null)
			pdfGraphsPane = new JTabbedPane();
		else
			while (pdfGraphsPane.getTabCount() > 0)
				pdfGraphsPane.removeTabAt(0);
		
		HistogramFunction[] aValExtras = null;
		HistogramFunction[] pValExtras = null;
		HistogramFunction[] logcValExtras = null;
		
		if (genericModel != null) {
			HistogramFunction genericA = genericModel.getPDF_a();
			HistogramFunction genericP = genericModel.getPDF_p();
			HistogramFunction genericLogC = genericModel.getPDF_logc();
			
			genericA.setName("Generic");
			genericP.setName("Generic");
			genericLogC.setName("Generic");

			
			HistogramFunction bayesianA;
			HistogramFunction bayesianP;
			HistogramFunction bayesianLogC;
			if (bayesianModel != null) {
				
				bayesianA = bayesianModel.getPDF_a();
				bayesianP = bayesianModel.getPDF_p();
				bayesianLogC = bayesianModel.getPDF_logc();
				
				if(bayesianA != null)
					bayesianA.setName("Bayesian");
				if(bayesianP != null)
					bayesianP.setName("Bayesian");
				if(bayesianLogC != null)
					bayesianLogC.setName("Bayesian");
				
				aValExtras = new HistogramFunction[] { genericA, bayesianA };
				pValExtras = new HistogramFunction[] { genericP, bayesianP };
				logcValExtras = new HistogramFunction[] { genericLogC, bayesianLogC };
			} else {
				aValExtras = new HistogramFunction[] { genericA };
				pValExtras = new HistogramFunction[] { genericP };
				logcValExtras = new HistogramFunction[] { genericLogC };
			}
		}
		
		HistogramFunction testPdf = seqSpecModel.getPDF_a();
		System.out.println("aHistMax = " + testPdf.getMaxY());
		testPdf = seqSpecModel.getPDF_p();
		System.out.println("pHistMax = " + testPdf.getMaxY());
		
		add1D_PDF(seqSpecModel.getPDF_a(), "a-value", aValExtras);
		add1D_PDF(seqSpecModel.getPDF_p(), "p-value", pValExtras);
		add1D_PDF(seqSpecModel.getPDF_logc(), "logc-value", logcValExtras);
		add2D_PDF(seqSpecModel.get2D_PDF_for_a_and_logc(), "a-value", "logc-value");
		add2D_PDF(seqSpecModel.get2D_PDF_for_a_and_p(), "a-value", "p-value");
		add2D_PDF(seqSpecModel.get2D_PDF_for_logc_and_p(), "logc-value", "p-value");
		
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
	
	private void plotExpectedAfershockMFDs(CalcProgressBar progress) {
		Double minDays = forecastStartTimeParam.getValue();
		validateParameter(minDays, "start time");
		Double maxDays = forecastEndTimeParam.getValue();
		validateParameter(maxDays, "end time");
		
		double minMag;
//		if (mainshock.getMag() < 6)
//			minMag = 3d;
//		else
//			minMag = 3d;
		
		minMag = Math.min(3d, mcParam.getValue());
		
		double maxMag = 9d;
		double deltaMag = 0.1;
		int numMag = (int)((maxMag - minMag)/deltaMag + 1.5);
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		List<ETAS_AftershockModel> models = Lists.newArrayList();
		List<String> names = Lists.newArrayList();
		List<Color> colors = Lists.newArrayList();
		
		if (seqSpecModel != null) {
			models.add(seqSpecModel);
			names.add("Seq. Specific");
			colors.add(sequence_specific_color);
		}
		if (genericModel != null) {
			models.add(genericModel);
			names.add("Generic");
			colors.add(generic_color);
		}
		if (bayesianModel != null) {
			models.add(bayesianModel);
			names.add("Bayesian");
			colors.add(bayesian_color);
		}
				
		double[] fractiles = { 0.025, 0.975 };
		
		for (int i=0; i<models.size(); i++) {
			String name = names.get(i);
			if (progress != null)
				progress.updateProgress(i, models.size(), "Calculating "+name+"...");
			
			System.out.println(name);
			
			ETAS_AftershockModel model = models.get(i);
			EvenlyDiscretizedFunc mode = model.getModalCumNumMFD(minMag, maxMag, numMag, minDays, maxDays);
			mode.setName(name+" Mode p50%");
//			EvenlyDiscretizedFunc mean = model.getMeanCumNumMFD(minMag, maxMag, numMag, minDays, maxDays);
//			mean.setName("Mean");
			Color c = colors.get(i);
			
			funcs.add(mode);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
			
//			funcs.add(mean);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, c));
			
//			EvenlyDiscretizedFunc[] fractilesFuncs = model.getCumNumMFD_FractileWithAleatoryVariability(
//					fractiles, minMag, maxMag, numMag, minDays, maxDays);
			
			
			for (int j= 0; j<fractiles.length; j++) {
			
				double f = fractiles[j];
				EvenlyDiscretizedFunc fractile = model.getCumNumMFD_Fractile(
						f, minMag, maxMag, numMag, minDays, maxDays);
				
//				EvenlyDiscretizedFunc fractile = fractilesFuncs[j];
				fractile.setName("p"+(float)(f*100d)+"%");
				
				funcs.add(fractile);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, c));
			}
		}
		
		if (progress != null) {
			progress.setIndeterminate(true);
			progress.setProgressMessage("Plotting...");
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
		
//		funcs.add(mainshockFunc);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
//		funcs.add(bathsFunc);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
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
			names.add("Bayesian");
		}

		if (seqSpecModel != null){
			models.add(seqSpecModel);
			names.add("Seq. Specific");
		}
		
		if (genericModel != null) {
			models.add(genericModel);
			names.add("Generic");
		}	
	
		
		GregorianCalendar eventDate = mainshock.getOriginTimeCal();
		GregorianCalendar startDate = new GregorianCalendar();
		Double minDays = forecastStartTimeParam.getValue();
		validateParameter(minDays, "start time");
		double startTime = eventDate.getTime().getTime() + minDays*ProbabilityModelsCalc.MILLISEC_PER_DAY;
		startDate.setTimeInMillis((long)startTime);
		
		for (int i=0; i<models.size(); i++) {
			Stopwatch watch = Stopwatch.createStarted();
			ETAS_AftershockModel model = models.get(i);
			String name = names.get(i);
			
			if (progress != null)
				progress.updateProgress(i, models.size(), "Calculating "+name+"...");
			
			ETAS_USGS_AftershockForecast forecast = new ETAS_USGS_AftershockForecast(model, eventDate, startDate);
			JTable jTable = new JTable(forecast.getTableModel());
			jTable.getTableHeader().setFont(jTable.getTableHeader().getFont().deriveFont(Font.BOLD));
			forecastTablePane.addTab(name, jTable);
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
	
	private void plotRateModel2D(CalcProgressBar progress){
		if (forecastMapPane == null)
			forecastMapPane = new JTabbedPane();
		else
			while (forecastMapPane.getTabCount() > 0)
				forecastMapPane.removeTabAt(0);

		double spacing = 0.025; 	// grid spacing in degrees
		double stressDrop = 3.0; 	//MPa
		double mainshockFitDuration = 0.5; //days
		
		if (progress != null)
			progress.updateProgress(0, 1, "Calculating Rate Model...");
		
		// TODO progress inside loop in this method
		GriddedGeoDataSet genericRateModel = genericModel.getRateModel2D(spacing, stressDrop, mainshockFitDuration);

		System.out.println("rateModel min max: " + + genericRateModel.getMinZ() + " " + genericRateModel.getMaxZ());
		
		
		if (progress != null)
			progress.updateProgress(1,1);
		
		CPT cpt;
		try {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(genericRateModel.getMinZ(), genericRateModel.getMaxZ());
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		XYZPlotSpec spec = new XYZPlotSpec(genericRateModel, cpt, "Spatial Forecast", "Longitude", "Latitude", "Epicentral Rate Density (km-2");
		
		XYZGraphPanel xyzGP = new XYZGraphPanel();
		forecastMapPane.addTab("Aftershock Rate", null, xyzGP);
		GriddedRegion reg = genericRateModel.getRegion();
		double xDelta = reg.getLonSpacing();
		double yDelta = reg.getLatSpacing();
		xyzGP.drawPlot(spec, false, false,
				new Range(reg.getMinGridLon()-0.5*xDelta, reg.getMaxGridLon()+0.5*xDelta),
				new Range(reg.getMinGridLat()-0.5*yDelta, reg.getMaxGridLat()+0.5*yDelta));
		
		if (tabbedPane.getTabCount() == forecast_map_tab_index)
			tabbedPane.addTab("Forecast Maps", null, forecastMapPane,
					"Forcast Maps (rate, intensity)");
		else
			Preconditions.checkState(tabbedPane.getTabCount() > forecast_map_tab_index, "Plots added out of order");
	}
	
	
	
	
	private class CalcStep {
		
		private String title;
		private String progressMessage;
		private Runnable run;
		boolean runInEDT;
		
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
			SwingUtilities.invokeLater(new Runnable() {
				
				@Override
				public void run() {
					if (progress.isVisible())
						progress.setVisible(false);
					progress.setModalityType(ModalityType.APPLICATION_MODAL);
					progress.setVisible(true);
				}
			});
			for (final CalcStep step : steps) {
				SwingUtilities.invokeLater(new Runnable() {
					
					@Override
					public void run() {
						progress.setTitle(step.title);
						progress.setProgressMessage(step.progressMessage);
						progress.pack();
					}
				});
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
				} else {
					try {
						step.run.run();
					} catch (Throwable e) {
						exception = e;
						break;
					}
				}
			}
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
		
		final CalcProgressBar progress = new CalcProgressBar(this, "", "", false);
		progress.setIndeterminate(true);
		
		if (param == eventIDParam || param == dataStartTimeParam || param == dataEndTimeParam
				|| param == regionEditParam) {
			setEnableParamsPostFetch(false);
			genericModel = null;
			bayesianModel = null;
			seqSpecModel = null;
			if(forecastTablePane != null)
				while (forecastTablePane.getTabCount() > 0)
					forecastTablePane.removeTabAt(0);
			if(pdfGraphsPane != null)
				while (pdfGraphsPane.getTabCount() > 0)
					pdfGraphsPane.removeTabAt(0);
			
		} else if (param == regionTypeParam || param == regionCenterTypeParam) {
			updateRegionParamList(regionTypeParam.getValue(), regionCenterTypeParam.getValue());
			setEnableParamsPostFetch(false);
		} else if (param == fetchButton) {
			setEnableParamsPostFetch(false);
			CalcStep fetchStep = new CalcStep("Fetching Events",
					"Contacting USGS ComCat Webservice. This is occasionally slow. "
					+ "If it fails, trying again often works.", new Runnable() {
						
						@Override
						public void run() {
							fetchEvents();
						}
					});
			CalcStep postFetchPlotStep = new CalcStep("Plotting Events/Data", "...", new Runnable() {
						
						@Override
						public void run() {
							doPostFetchPlots();
						}
					});
			CalcRunnable run = new CalcRunnable(progress, fetchStep, postFetchPlotStep);
			new Thread(run).start();
			
			
			
		} else if (param == loadCatalogButton) {
			if (loadCatalogChooser == null)
				loadCatalogChooser = new JFileChooser();
			int ret = loadCatalogChooser.showOpenDialog(this);
			if (ret == JFileChooser.APPROVE_OPTION) {
				final File file = loadCatalogChooser.getSelectedFile();
				CalcStep loadStep = new CalcStep("Loading events", "...", new Runnable() {
					
					@Override
					public void run() {
						try {
							loadCatalog(file);
						} catch (IOException e) {
							ExceptionUtils.throwAsRuntimeException(e);
						}
					}
				});
				CalcStep postFetchPlotStep = new CalcStep("Plotting Events/Data", "...", new Runnable() {
					
					@Override
					public void run() {
						doPostFetchPlots();
					}
				});
				CalcRunnable run = new CalcRunnable(progress, loadStep, postFetchPlotStep);
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
			
			// if a generic model has already been computed, update it to reflect new Mc/Mref productivity adjustment
			if(genericModel != null){
				System.out.println("Adjusting generic productivity to reflect Mc change...");
				genericModel = new ETAS_AftershockModel_Generic(mainshock, aftershocks, genericParams, 
					dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
					forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), mcParam.getValue(),
					9.5, 100, 0);
			}
			
				
			
		} else if (param == computeBButton) {
			setEnableParamsPostComputeB(false);
			CalcStep bStep = new CalcStep("Computing b", "...", new Runnable() {
				
				@Override
				public void run() {
					Double mc = mcParam.getValue();
					validateParameter(mc, "Mc");
					
					double magPrecision = magPrecisionParam.getValue();
					
					ObsEqkRupList filteredRupList = aftershocks.getRupsAboveMag(mc);
					double b = ETAS_StatsCalc.getMaxLikelihood_b_value(filteredRupList, mc, magPrecision);
					System.out.println("Num rups ≥ Mc = "+filteredRupList.size());
					System.out.println("Computed b-value: "+b);
					bParam.setValue(b);
					bParam.getEditor().refreshParamEditor();
					
					setEnableParamsPostComputeB(true);
					
					tabbedPane.setSelectedIndex(mag_num_tab_index);
					
					
				}
			});
			CalcRunnable run = new CalcRunnable(progress, bStep);
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
//		} else if (param == gParam || param == hParam || param == mCatParam) {
		} else if (param == rmaxParam){
			setEnableParamsPostAfershockParams(false);
		} else if (param == computeAftershockParamsButton) {
			setEnableParamsPostAfershockParams(false);
			CalcStep computeStep = new CalcStep("Computing Aftershock Params", "...", new Runnable() {

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
//						Double g = gParam.getValue();
//						validateParameter(g, "G");
//						
//						Double h = hParam.getValue();
//						validateParameter(h, "H");
//						
//						Double mCat = mCatParam.getValue();
//						validateParameter(mCat, "Mcat");
//						
						Double rmax = rmaxParam.getValue();
						validateParameter(rmax, "rmax");
//						
//						model = new RJ_AftershockModel_SequenceSpecific(mainshock, aftershocks, mCat, g, h, b,
//								dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
//								aRange.getLowerBound(), aRange.getUpperBound(), aNum,
//								pRange.getLowerBound(), pRange.getUpperBound(), pNum,
//								cRange.getLowerBound(), cRange.getUpperBound(), cNum);
						
						seqSpecModel = new ETAS_AftershockModel_SequenceSpecific(mainshock, aftershocks, rmax,
								ETAS_StatsCalc.linspace(aRange.getLowerBound(), aRange.getUpperBound(), aNum),
								ETAS_StatsCalc.linspace(pRange.getLowerBound(), pRange.getUpperBound(), pNum),
								ETAS_StatsCalc.logspace(cRange.getLowerBound(), cRange.getUpperBound(), cNum),
								1, b, magRefParam.getValue(),	dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
								forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(),	// forecast start time, forecast end time 
								9.5, 100, 0);	//max sim mag, max number of generations, number of simulations
								
						
						
					} else {
						seqSpecModel = new ETAS_AftershockModel_SequenceSpecific(mainshock, aftershocks,
								ETAS_StatsCalc.linspace(aRange.getLowerBound(), aRange.getUpperBound(), aNum),
								ETAS_StatsCalc.linspace(pRange.getLowerBound(), pRange.getUpperBound(), pNum),
								ETAS_StatsCalc.logspace(cRange.getLowerBound(), cRange.getUpperBound(), cNum),
								1, b, magRefParam.getValue(), dataStartTimeParam.getValue(), dataEndTimeParam.getValue(), 
								forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(),
									9.5, 100, 0);
								
					}
				}
				
			});
			CalcStep plotStep = new CalcStep("Plotting Model PDFs", "...", new Runnable() {
				
				@Override
				public void run() {
					aValParam.setValue(seqSpecModel.getMaxLikelihood_a());
					aValParam.getEditor().refreshParamEditor();
					pValParam.setValue(seqSpecModel.getMaxLikelihood_p());
					pValParam.getEditor().refreshParamEditor();
					cValParam.setValue(seqSpecModel.getMaxLikelihood_c());
					cValParam.getEditor().refreshParamEditor();
					
					bayesianModel = null;
					if (genericModel != null) {
						bayesianModel = new ETAS_AftershockModel_Bayesian(seqSpecModel, genericModel);
					}
										
					plotPDFs();
					setEnableParamsPostAfershockParams(true);
					plotCumulativeNum();
					tabbedPane.setSelectedIndex(pdf_tab_index);
				}
			}, true);
			CalcRunnable run = new CalcRunnable(progress, computeStep, plotStep);
			new Thread(run).start();
		} else if (param == computeAftershockForecastButton) {
			CalcStep plotStep = new CalcStep("Computing/Plotting Forecast MFDs", "This can take some time...",
					new Runnable() {

						@Override
						public void run() {
							System.out.println("Generic model");
							genericModel.setMagComplete(mcParam.getValue());
							genericModel.computeNewForecast(dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
									forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), 10000);
							
							seqSpecModel.setMagComplete(mcParam.getValue());
							System.out.println("Sequence Specific model");
							seqSpecModel.computeNewForecast(dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
									forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), 10000);
							
							bayesianModel.setMagComplete(mcParam.getValue());
							System.out.println("Bayesian model");
							bayesianModel.computeNewForecast(dataStartTimeParam.getValue(), dataEndTimeParam.getValue(),
									forecastStartTimeParam.getValue(), forecastEndTimeParam.getValue(), 10000);
							
							
							Stopwatch watch = Stopwatch.createStarted();
							plotExpectedAfershockMFDs(progress);
							watch.stop();
							System.out.println(
									"Took "+watch.elapsed(TimeUnit.SECONDS)+"s to compute/plot aftershock MFDs");
							
							//update cumNumPlot
							plotCumulativeNum();
							
						}
				
			});
			CalcStep tableStep = new CalcStep("Computing Forecast Table", "This can take some time...",
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
			});	
			CalcStep mapStep = new CalcStep("Computing Rate Map", "This can take some time...",
					new Runnable() {

						@Override
						public void run() {
							Stopwatch watch = Stopwatch.createStarted();

							plotRateModel2D(progress);
							
							watch.stop();
							System.out.println(
									"Took "+watch.elapsed(TimeUnit.SECONDS)+"s to compute/plot rate map");
							//								tabbedPane.setSelectedIndex(aftershock_expected_index);
						}
			});
			CalcRunnable run = new CalcRunnable(progress, plotStep, tableStep, mapStep);
			new Thread(run).start();
		
		}
	}
	
	private void doPostFetchCalculations() {
		aftershockMND = ObsEqkRupListCalc.getMagNumDist(aftershocks, 1.05, 81, 0.1);
		mmaxc = ETAS_StatsCalc.getMmaxC(aftershockMND);
		mcParam.setValue(mmaxc+0.5);
		mcParam.getEditor().refreshParamEditor();
		
//		//update data end time with the maximum time in the catalog
//		ObsEqkRupture lastAftershock = aftershocks.get(aftershocks.size()-1);
//		System.out.println(lastAftershock.getOriginTime());
//		System.out.println(mainshock.getOriginTime());
//		dataEndTimeParam.setValue((double) (lastAftershock.getOriginTime() - mainshock.getOriginTime()) / (double) ETAS_StatsCalc.MILLISEC_PER_DAY );
//		dataEndTimeParam.getEditor().refreshParamEditor();
//		
		
	}
	
	private void doPostFetchPlots() {
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
		
		rmaxParam.getEditor().setEnabled(enabled && timeDepMcParam.getValue());
//		gParam.getEditor().setEnabled(enabled && timeDepMcParam.getValue());
//		hParam.getEditor().setEnabled(enabled && timeDepMcParam.getValue());
//		mCatParam.getEditor().setEnabled(enabled && timeDepMcParam.getValue());
		
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
		forecastStartTimeParam.getEditor().setEnabled(enabled);
		forecastEndTimeParam.getEditor().setEnabled(enabled);
		computeAftershockForecastButton.getEditor().setEnabled(enabled);
		if (!enabled)
			seqSpecModel = null;
	}
	
	public static void main(String[] args) {
		new AftershockStatsGUI_ETAS();
	}

}
