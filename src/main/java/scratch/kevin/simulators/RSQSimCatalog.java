package scratch.kevin.simulators;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.dom4j.Attribute;
import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.metadata.XMLSaveable;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.XMLUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RupSetDeformationModel;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.imr.attenRelImpl.ngaw2.FaultStyle;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.distCalc.SimRuptureDistCalcUtils.LocationElementDistanceCacheFactory;
import org.opensha.sha.simulators.iden.CatalogLengthLoadIden;
import org.opensha.sha.simulators.iden.EventIDsRangeIden;
import org.opensha.sha.simulators.iden.EventIDsRupIden;
import org.opensha.sha.simulators.iden.EventTimeIdentifier;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.iden.SectionIDIden;
import org.opensha.sha.simulators.iden.SkipYearsLoadIden;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader.TransVersion;
import org.opensha.sha.simulators.srf.RSQSimTransValidIden;
import org.opensha.sha.simulators.utils.RSQSimEqkRupture;
import org.opensha.sha.simulators.utils.RSQSimSubSectEqkRupture;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SlipAlongSectAlgorithm;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simulators.nz.NZ_CompModels.NZ_CompDefModel;
import scratch.kevin.simulators.nz.NZ_CompModels.NZ_CompFaultModel;
import scratch.kevin.simulators.plots.AbstractPlot;
import scratch.kevin.simulators.plots.MFDPlot;
import scratch.kevin.simulators.plots.MagAreaScalingPlot;
import scratch.kevin.simulators.plots.MomentRateVaribilityPlot;
import scratch.kevin.simulators.plots.NormalizedFaultRecurrenceIntervalPlot;
import scratch.kevin.simulators.plots.PaleoOpenIntervalPlot;
import scratch.kevin.simulators.plots.PaleoRecurrencePlot;
import scratch.kevin.simulators.plots.RecurrenceIntervalPlot;
import scratch.kevin.simulators.plots.RuptureVelocityPlot;
import scratch.kevin.simulators.plots.SectParticipationNucleationPlot;
import scratch.kevin.simulators.plots.SectionRecurrenceComparePlot;
import scratch.kevin.simulators.plots.SectionRecurrenceComparePlot.SectType;
import scratch.kevin.simulators.plots.SlipAlongRupturePlot;
import scratch.kevin.simulators.plots.SlipLengthScalingPlot;
import scratch.kevin.simulators.plots.SlipRateComparePlot;
import scratch.kevin.simulators.plots.StationarityPlot;
import scratch.kevin.simulators.plots.ElasticReboundTriggeringPlot;
import scratch.kevin.simulators.plots.U3StyleNormalizedRuptureRecurrenceIntervalPlot;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.FilterMethod;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.rotation.RSQSimRotatedRupVariabilityMagDistPageGen.RuptureType;

public class RSQSimCatalog implements XMLSaveable {
	
	public enum Catalogs {
		BRUCE_2142("bruce/rundir2142", "Bruce 2142", "Bruce Shaw", cal(2017, 6, 16),
				"Old projection; slip weakening; stress loaded; no creep correction",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2194("bruce/rundir2194", "Bruce 2194", "Bruce Shaw", cal(2017, 7, 5),
				"Catalog with decent large event scaling and distribution of sizes while not using"
				+ " any of the enhanced frictional weakening terms.", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_UCERF3_millionElement("JG_UCERF3_millionElement", "U3 1mil Element Test", "Jacqui Gilchrist", cal(2017, 9, 27),
				"Test 1 million element catalog on UCERF3 fault system, ~0.25 km^2 trianglar elements",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2194_LONG("rundir2194_long", "Bruce 2194 Long", "Bruce Shaw (extended by Jacqui Gilchrist)", cal(2017, 8, 31),
				"Catalog with decent large event scaling and distribution of sizes while not using"
				+ " any of the enhanced frictional weakening terms.", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2273("bruce/rundir2273", "Bruce 2273", "Bruce Shaw", cal(2017, 10, 13),
				"Stress loading, more refined geometry, does not contain projection fix (some location discrepancies "
				+ "are present relative to UCERF3 faults).", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2310("bruce/rundir2310", "Bruce 2310", "Bruce Shaw", cal(2017, 10, 16),
				"Backslip loading, more refined geometry, projection fix (but all faults surface breaking)",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2320("bruce/rundir2320", "Bruce 2320", "Bruce Shaw", cal(2017, 10, 17),
				"Backslip loading, less refined geometry, projection fix (but all\n" + 
				"faults surface breaking), same as rundir2310 but less resolved",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2336("bruce/rundir2336", "Bruce 2336", "Bruce Shaw", cal(2017, 10, 20),
				"Larger slip velocity (1.5 m/s), backslipFromStress loading",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2337("bruce/rundir2337", "Bruce 2337", "Bruce Shaw", cal(2017, 10, 20),
				"Larger slip velocity (2.0 m/s), backslipFromStress loading",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2326("bruce/rundir2326", "Bruce 2326", "Bruce Shaw", cal(2017, 10, 23),
				"reference_1: a=.001 b=.008  Veq=1.0  sigmaN=100. backslipFromStress",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2342("bruce/rundir2342", "Bruce 2342", "Bruce Shaw", cal(2017, 10, 23),
				"larger Veq=1.2          relative to reference_1",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2343("bruce/rundir2343", "Bruce 2343", "Bruce Shaw", cal(2017, 10, 23),
				"smaller Veq=0.8        relative to reference_1",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2349("bruce/rundir2349", "Bruce 2349", "Bruce Shaw", cal(2017, 10, 23),
				"smaller sigmaN=80   relative to reference_1",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_2194_K2("rundir2194_K2", "JG 2194 K2", "Jacqui Gilchrist", cal(2017, 10, 16),
				"Keith's fault geometry, normal backslip with U3 geologic long-term slip rates,"
				+ " and the same parameter values as Bruce's 2194",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_tuneBase1m("tuneBase1m", "JG Tune Base 1M", "Jacqui Gilchrist", cal(2017, 11, 2),
				"U3 fault geometry with 1km^2 triangles, normal backslip loading with U3 geologic slip rates,"
				+ "calibrated to U3 supraseismogenic recurrence intervals, and default a/b",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_modLoad_testB("modLoad_testB", "JG Mod Load Test B", "Jacqui Gilchrist", cal(2017, 11, 14),
				"Bruce's modified loading with higher values of frictional parameters",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_baseCatalog2("baseCatalog2", "JG Base Catalog 2", "Jacqui Gilchrist", cal(2017, 11, 16),
				"Untuned version of tuneBase1m. Same fault model and frictional parameters, without any stress adjustments",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_tuneBaseSW_1e5("tuneBaseSW_1e5", "JG Tune Base SW", "Jacqui Gilchrist", cal(2017, 11, 20),
				"Tuned, additional slip weakening parameters using Keith's fault geometry. muSlipAmp = 0.2, muSlipInvDist_1 = 2.0, cohesion = 6.",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_baseCatalogSW_10("baseCatalogSW_10", "JG Base SW", "Jacqui Gilchrist", cal(2017, 11, 20),
				"Untuned, additional slip weakening parameters using Keith's fault geometry. muSlipAmp = 0.2, muSlipInvDist_1 = 2.0, cohesion = 6.",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_tunedBase1m_ddotEQmod("tunedBase1m_ddotEQmod", "JG Tune Base Mod Vel", "Jacqui Gilchrist", cal(2017, 11, 26),
				"New version of tuneBase1m, with patch-specific slip velocities.",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2381("bruce/rundir2381", "Bruce 2381", "Bruce Shaw", cal(2017, 12, 22),
				"fracCreep=0.5", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2388("bruce/rundir2388", "Bruce 2388", "Bruce Shaw", cal(2017, 12, 22),
				"fracCreep=0.25", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2410("bruce/rundir2410", "Bruce 2410", "Bruce Shaw", cal(2017, 12, 22),
				"factorNormal=2.0;  50ky reference time; 212ky spinup", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2411("bruce/rundir2411", "Bruce 2411", "Bruce Shaw", cal(2017, 12, 22),
				"factorNormal=2.0;  50ky reference time; 85ky spinup", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2412("bruce/rundir2412", "Bruce 2412", "Bruce Shaw", cal(2017, 12, 22),
				"factorNormal=2.0;  50ky reference time; 85ky spinup; srt(slipRate)", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2413("bruce/rundir2413", "Bruce 2413", "Bruce Shaw", cal(2017, 12, 22),
				"factorNormal=2.0;  100ky reference time; 85ky spinup; sqrt(slipRate)", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2457("bruce/rundir2457", "Bruce 2457", "Bruce Shaw", cal(2018, 1, 14),
				"new loading;  fCreep=0.25", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2495("bruce/rundir2495", "Bruce 2495", "Bruce Shaw", cal(2018, 1, 29),
				"flat loaded.  fracCreep=0.5.  maxDepth=14.4", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2579("bruce/rundir2579", "Bruce 2579", "Bruce Shaw", cal(2018, 2, 07),
				"straight loaded;  fracCreep=0.5;  H=18 (2,12,4);  stressMult=1.2;  neighbors",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2585("bruce/rundir2585", "Bruce 2585", "Bruce Shaw", cal(2018, 2, 10),
				"Longer run; else same as r2579. straight loaded;  fracCreep=0.5; H=18 (2,12,4); stressMult=1.2; neighbors",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2592("bruce/rundir2592", "Bruce 2592", "Bruce Shaw", cal(2018, 2, 11),
				"straight loaded;  fracCreep=0.5;  H=16 (2,11,3); stressMult=1.2; neighbors",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2620("bruce/rundir2620", "Bruce 2620", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   b=.006", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2621("bruce/rundir2621", "Bruce 2621", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   b=.007", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2622("bruce/rundir2622", "Bruce 2622", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   b=.009", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2623("bruce/rundir2623", "Bruce 2623", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   b=.010", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2624("bruce/rundir2624", "Bruce 2624", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   a=.0008", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2625("bruce/rundir2625", "Bruce 2625", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   a=.0009", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2626("bruce/rundir2626", "Bruce 2626", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   a=.0011", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2627("bruce/rundir2627", "Bruce 2627", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   a=.0012", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2628("bruce/rundir2628", "Bruce 2628", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   N=80", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2629("bruce/rundir2629", "Bruce 2629", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   N=90", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2630("bruce/rundir2630", "Bruce 2630", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   N=110", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2631("bruce/rundir2631", "Bruce 2631", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   N=120", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2585_1MYR("rundir2585_1myr", "Bruce 2585 1myr", "Bruce Shaw/Jacqui Gilchrist", cal(2018, 3, 12),
				"Extended version of Bruce's 2585 to 1 million years",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2616("bruce/rundir2616", "Bruce 2616", "Bruce Shaw", cal(2018, 3, 19),
				"similar to r2585, but bigger seismogenic depth: H=18 (2,13,3)",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2636("bruce/rundir2636", "Bruce 2636", "Bruce Shaw", cal(2018, 3, 28),
				"sensitivity test, diff r2585 a=.0013",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2637("bruce/rundir2637", "Bruce 2637", "Bruce Shaw", cal(2018, 3, 28),
				"sensitivity test, diff r2585  N=130",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2665("bruce/rundir2665", "Bruce 2665", "Bruce Shaw", cal(2018, 4, 16),
				"dx/2, LatCut=37, rateCut=0.2mm/yr, interpolated nearest",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2666("bruce/rundir2666", "Bruce 2666", "Bruce Shaw", cal(2018, 4, 17),
				"dx/4, LatCut=37, rateCut=2.0mm/yr, interpolated nearest",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2667("bruce/rundir2667", "Bruce 2667", "Bruce Shaw", cal(2018, 4, 23),
				"dx/4, LatCut=37, rateCut=2.0mm/yr, interpolated nearest, like r2666 but larger b=.015",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2734("bruce/rundir2734", "Bruce 2734", "Bruce Shaw", cal(2018, 6, 5),
				"Finite receiver patch fracArea=0.8, all else same as r2585",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2737("bruce/rundir2737", "Bruce 2737", "Bruce Shaw", cal(2018, 6, 5),
				"Finite receiver patch fracArea=0.9, all else same as r2585",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2740("bruce/rundir2740", "Bruce 2740", "Bruce Shaw", cal(2018, 6, 27),
				"FinitePatch fracArea=0.8, 48Hr, all else same as r2585",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2742("bruce/rundir2742", "Bruce 2742", "Bruce Shaw", cal(2018, 6, 27),
				"FinitePatch fracArea=0.85,  all else same as r2585",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2743("bruce/rundir2743", "Bruce 2743", "Bruce Shaw", cal(2018, 6, 27),
				"FinitePatch fracArea=0.8,  b=.007",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2744("bruce/rundir2744", "Bruce 2744", "Bruce Shaw", cal(2018, 6, 27),
				"FinitePatch fracArea=0.8,  b=.009",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2768("bruce/rundir2768", "Bruce 2768", "Bruce Shaw", cal(2018, 8, 14),
				"FinitePatch fracArea=0.90, 48Hr, all else same as r2585",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2769("bruce/rundir2769", "Bruce 2769", "Bruce Shaw", cal(2018, 8, 14),
				"FinitePatch fracArea=0.95, 48Hr, all else same as r2585",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
//		JG_SWEEPTEST_2763V2("test25a589d_2763V2", "JG Test Rerun 2763", "Jacqui Gilchrist", cal(2018, 8, 21),
//				"Rerun test of Bruce's rundir2763, reading new stiffness matrices",
//				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
//		JG_SWEEPTEST_J("test25a589d_J", "JG Test Calc Siffness", "Jacqui Gilchrist", cal(2018, 8, 21),
//				"Test of calculating and writing new stiffness matrices",
//				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
//		JG_SWEEPTEST_J2("test25a589d_J2", "JG Test Calc Siffness 2", "Jacqui Gilchrist", cal(2018, 8, 22),
//				"Test of reading in the stiffness files written by previous run test25a589d_J using my executable",
//				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
//		JG_SWEEPTEST_2763V3("test25a589d_flt2763V3", "JG Test Rerun 2763 V3", "Jacqui Gilchrist", cal(2018, 8, 22),
//				"Test of reading in the stiffness files written by previous run test25a589d_J using Bruce's executable",
//				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
//		JG_SWEEPTEST_MASTER_EXEC("defaultModel_JfM", "JG Test Master Exec", "Jacqui Gilchrist", cal(2018, 8, 23),
//				"Calculating stiffness files using the executable from the master branch",
//				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2829("bruce/rundir2829", "Bruce 2829", "Bruce Shaw", cal(2018, 10, 18),
				"fracArea=0.95 ; NEW variableSpeed  s2ddf=.9 ddfmin=.1;  b=.01 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_3013("bruce/rundir3013", "Bruce 3013", "Bruce Shaw", cal(2019, 1, 13),
				"New const dip fault smoothing.  smoothF=1e6.  V=1.  fracArea=0.99.  b=.01",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_3014("bruce/rundir3014", "Bruce 3014", "Bruce Shaw", cal(2019, 1, 13),
				"New const dip fault smoothing.  smoothF=1e6.  V=1.  fracArea=0.99.  b=.013",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_3032("bruce/rundir3032", "Bruce 3032", "Bruce Shaw", cal(2019, 1, 19),
				"SmoothF=1e5.   Unconnected discontinuities. V=1.  fracArea=0.99. b=.011",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_3062("bruce/rundir3062", "Bruce 3062", "Bruce Shaw", cal(2019, 2, 8),
				"SmoothF=1e5. Connected discontinuities. V=1. fracArea=0.99. b=.011",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_3065("bruce/rundir3065", "Bruce 3065", "Bruce Shaw", cal(2019, 2, 8),
				"Variable normal stress dsigdsA =0.6.  Rest same as r3062",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_3067("bruce/rundir3067", "Bruce 3067", "Bruce Shaw", cal(2019, 2, 10),
				"Variable normal stress dsigdsA =0.666.  Rest same as r3062",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_3164("bruce/rundir3164", "Bruce 3164", "Bruce Shaw", cal(2019, 3, 9),
				"AReduceDelay  tCausalFactor=3.0 areaFrac=0.99  V=1  b=.027 a=.015 fA=.005 ; Smooth Model",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_3165("bruce/rundir3165", "Bruce 3165", "Bruce Shaw", cal(2019, 3, 18),
				"AReduceDelay  tCausalFactor=3.0 areaFrac=0.99  V=1  b=.017 a=.005 fA=.005 ; Smooth Model",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_3271("bruce/rundir3271", "Bruce 3271", "Bruce Shaw", cal(2019, 5, 29),
				"a=.006 b=.018  more b=1 like",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		USEIT_600K_MOD_PARAMS("useit_longA_600k", "UseIT 600k Mod Params", "UseIT FAST Interns", cal(2019, 7, 7),
				"a=.006 b=.018, all else 2585",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4312("bruce/rundir4312", "Bruce 4312", "Bruce Shaw", cal(2019, 7, 24),
				"same as r2585 but variable slip speed.  fracArea=0 ; variableSpeed s2ddf=.9 ddfmin=.1; b=.008 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4313("bruce/rundir4313", "Bruce 4313", "Bruce Shaw", cal(2019, 7, 25),
				"same as r2585 but variable slip speed.  fracArea=0 ; variableSpeed s2ddf=.9 ddfmin=.1; b=.007 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4314("bruce/rundir4314", "Bruce 4314", "Bruce Shaw", cal(2019, 7, 25),
				"same as r2585 but variable slip speed.  fracArea=0 ; variableSpeed s2ddf=.9 ddfmin=.1; b=.006 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4316("bruce/rundir4316", "Bruce 4316", "Bruce Shaw", cal(2019, 7, 26),
				"same as r2585 but variable slip speed.  fracArea=0 ; variableSpeed s2ddf=.8 ddfmin=.1; b=.006 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4317("bruce/rundir4317", "Bruce 4317", "Bruce Shaw", cal(2019, 7, 26),
				"variable slip speed.  fracArea=0.8 ; variableSpeed s2ddf=.8 ddfmin=.1; b=.008 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4320("bruce/rundir4320", "Bruce 4320", "Bruce Shaw", cal(2019, 7, 28),
				"variable slip speed. fracArea=0.7 ; variableSpeed s2ddf=.7 ddfmin=.01; b=.008 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4322("bruce/rundir4322", "Bruce 4322", "Bruce Shaw", cal(2019, 7, 28),
				"variable slip speed. fracArea=0.8 ; variableSpeed s2ddf=.7 ddfmin=.01; b=.010 a=.003",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4576("bruce/rundir4576", "Bruce 4576", "Bruce Shaw", cal(2019, 10, 16),
				"same as r2585 but larger distanceTop=4.0; H=18 (4,10,4)",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4655("bruce/rundir4655", "Bruce 4655", "Bruce Shaw", cal(2019, 12, 9),
				"variable slip speed. fracArea=0 ; variableSpeed s2ddf=.8 ddfmin=.1; b=.008 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4660("bruce/rundir4660", "Bruce 4660", "Bruce Shaw", cal(2019, 12, 11),
				"variable slip speed. fracArea=0 ; variableSpeed s2ddf=.8 ddfmin=.1; b=.007 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4661("bruce/rundir4661", "Bruce 4661", "Bruce Shaw", cal(2019, 12, 11),
				"variable slip speed. fracArea=0 ; variableSpeed s2ddf=.8 ddfmin=.1; b=.006 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4682("bruce/rundir4682", "Bruce 4682", "Bruce Shaw", cal(2019, 12, 12),
				"higher resolution delta=1.25; variable slip speed. fracArea=0 ; "
				+ "variableSpeed s2ddf=.8 ddfmin=.1; b=.008 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4665("bruce/rundir4665", "Bruce 4665", "Bruce Shaw", cal(2019, 12, 18),
				"Long run. variable slip speed. fracArea=0 ; variableSpeed s2ddf=.8 ddfmin=.1; b=.007 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4667("bruce/rundir4667", "Bruce 4667", "Bruce Shaw", cal(2019, 12, 23),
				"Long run. variable slip speed. fracArea=0 ; variableSpeed s2ddf=.8 ddfmin=.1; b=.007 a=.001;"
				+ " higher tau0 to reduce transient",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4780("bruce/rundir4780", "Bruce 4780", "Bruce Shaw", cal(2020, 1, 23),
				"draft causal.  tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=.4 ; b=.01 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4812("bruce/rundir4812", "Bruce 4812", "Bruce Shaw", cal(2020, 1, 24),
				"draft causal.  tCausalF=0.65 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=.4 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4819("bruce/rundir4819", "Bruce 4819", "Bruce Shaw", cal(2020, 1, 29),
				"draft causal.  tCausalF=0.65 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=.4 ; b=.009 a=.001"
				+ " . [V>0 constraint]",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4820("bruce/rundir4820", "Bruce 4820", "Bruce Shaw", cal(2020, 1, 29),
				"variable slip speed. fracArea=0 ; variableSpeed s2ddf=.8 ddfmin=.1; b=.007 a=.001"
				+ " . [V>0 constraint]",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4827("bruce/rundir4827", "Bruce 4827", "Bruce Shaw", cal(2020, 1, 30),
				"draft causal.  tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=.4 ; b=.0085 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4836("bruce/rundir4836", "Bruce 4836", "Bruce Shaw", cal(2020, 2, 3),
				"draft causal.  tCausalF=0.65 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=.4 ; b=.0085 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4838("bruce/rundir4838", "Bruce 4838", "Bruce Shaw", cal(2020, 2, 3),
				"draft causal.  tCausalF=0.74 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=.4 ; b=.0085 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4841("bruce/rundir4841", "Bruce 4841", "Bruce Shaw", cal(2020, 2, 3),
				"draft causal.  tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=.3 ; b=.0085 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
//		BRUCE_4851("bruce/rundir4851", "Bruce 4851", "Bruce Shaw", cal(2020, 2, 3),
//				"higher res.  draft causal.  tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=.4 "
//				+ "; b=.0085 a=.001",
//				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4853("bruce/rundir4853", "Bruce 4853", "Bruce Shaw", cal(2020, 2, 3),
				"higher res.  draft causal.  tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=.3 "
				+ "; b=.0085 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4848("bruce/rundir4848", "Bruce 4848", "Bruce Shaw", cal(2020, 2, 6),
				"draft causal.  tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.22 ; b=.0085 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4849("bruce/rundir4849", "Bruce 4849", "Bruce Shaw", cal(2020, 2, 6),
				"draft causal.  tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.1 ; b=.0085 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4860("bruce/rundir4860", "Bruce 4860", "Bruce Shaw", cal(2020, 2, 6),
				"higher res.  draft causal.  tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.3 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4862("bruce/rundir4862", "Bruce 4862", "Bruce Shaw", cal(2020, 2, 6),
				"higher res.  draft causal.  tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.1 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		TEST_DOUBLE_4860("test_double_4860", "Test Double 4860", "Kevin Milner", cal(2020, 2, 7),
				"This is a test of catalog and transition file combining code, containing two copies of 4860. The "
				+ "first 5000 years are also skipped from the second copy.",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4860_10X("rundir4860_multi_combine", "Bruce 4860 10X Combine", "Kevin Milner", cal(2020, 2, 12),
				"Stitch of 10 different seeds of 4860 (should not be used for time-dependent calculations)",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4930("bruce/rundir4930", "Bruce 4930", "Bruce Shaw", cal(2020, 2, 27),
				"higher res. draft causal. tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.25 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4935("bruce/rundir4935", "Bruce 4935", "Bruce Shaw", cal(2020, 2, 27),
				"higher res. draft causal. tCausalF=0.65 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.15 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4950("bruce/rundir4950", "Bruce 4950", "Bruce Shaw", cal(2020, 3, 31),
				"same as 4860, but with an I/O fix to (hopefully) correct the double peak in patch velocities",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4962("bruce/rundir4962", "Bruce 4962", "Bruce Shaw", cal(2020, 4, 2),
				"tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.3 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4963("bruce/rundir4963", "Bruce 4963", "Bruce Shaw", cal(2020, 4, 2),
				"tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.3 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4964("bruce/rundir4964", "Bruce 4964", "Bruce Shaw", cal(2020, 4, 2),
				"tCausalF=0.67 ; fracArea=0.9 ; varV s2ddf=.8 ddfmin=0.3 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4965("bruce/rundir4965", "Bruce 4965", "Bruce Shaw", cal(2020, 4, 2),
				"tCausalF=0.67 ; fracArea=0.95 ; varV s2ddf=.8 ddfmin=0.3 ; b=.010 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4970("bruce/rundir4970", "Bruce 4970", "Bruce Shaw", cal(2020, 4, 5),
				"tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.3 ; b=.010 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4976("bruce/rundir4976", "Bruce 4976", "Bruce Shaw", cal(2020, 4, 5),
				"tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.9 ddfmin=0.3 ; b=.010 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4977("bruce/rundir4977", "Bruce 4977", "Bruce Shaw", cal(2020, 4, 11),
				"tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.3 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4978("bruce/rundir4978", "Bruce 4978", "Bruce Shaw", cal(2020, 4, 11),
				"higher res.  tCausalF=0.7 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.3 ; b=.010 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4979("bruce/rundir4979", "Bruce 4979", "Bruce Shaw", cal(2020, 4, 11),
				"higher res.  tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.3 ; b=.010 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4983("bruce/rundir4983", "Bruce 4983", "Bruce Shaw", cal(2020, 4, 20),
				"tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.3 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4984("bruce/rundir4984", "Bruce 4984", "Bruce Shaw", cal(2020, 4, 20),
				"higher res.  tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.3 ; b=.010 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4983_STITCHED("rundir4983_stitched", "Bruce 4983 (Stitched)", "Bruce Shaw", cal(2020, 4, 27),
				"Stitched extension of 4983, 4 runs",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_4987("bruce/rundir4987", "Bruce 4987", "Bruce Shaw", cal(2020, 4, 29),
				"VeqMax=4.0; tCausalF=0.63 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.3 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5040("bruce/rundir5040", "Bruce 5040", "Bruce Shaw", cal(2020, 9, 17),
				"updated tdelay.  tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.3 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5041("bruce/rundir5041", "Bruce 5041", "Bruce Shaw", cal(2020, 9, 17),
				"updated tdelay.  tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.2 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5042("bruce/rundir5042", "Bruce 5042", "Bruce Shaw", cal(2020, 9, 17),
				"updated tdelay.  tCausalF=0.70 ; fracArea=0.9 ; varV s2ddf=.8 ddfmin=0.2 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5043("bruce/rundir5043", "Bruce 5043", "Bruce Shaw", cal(2020, 9, 17),
				"updated tdelay.  tCausalF=0.70 ; fracArea=0.9 ; varV s2ddf=.8 ddfmin=0.2 ; b=.010 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5044("bruce/rundir5044", "Bruce 5044", "Bruce Shaw", cal(2020, 9, 17),
				"updated tdelay.  tCausalF=0.70 ; fracArea=0.9 ; varV s2ddf=.8 ddfmin=0.1 ; b=.010 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5045("bruce/rundir5045", "Bruce 5045", "Bruce Shaw", cal(2020, 9, 28),
				"updated tdelay.  tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.1  ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5046("bruce/rundir5046", "Bruce 5046", "Bruce Shaw", cal(2020, 9, 28),
				"updated tdelay.  tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.25 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5047("bruce/rundir5047", "Bruce 5047", "Bruce Shaw", cal(2020, 9, 28),
				"updated tdelay.  tCausalF=0.67 ; fracArea=0.8 ; varV s2ddf=.8 ddfmin=0.15 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5048("bruce/rundir5048", "Bruce 5048", "Bruce Shaw", cal(2020, 9, 28),
				"updated tdelay.  tCausalF=0.67 ; fracArea=0.85 ; varV s2ddf=.8 ddfmin=0.2 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5049("bruce/rundir5049", "Bruce 5049", "Bruce Shaw", cal(2020, 9, 28),
				"updated tdelay.  tCausalF=0.67 ; fracArea=0.75 ; varV s2ddf=.8 ddfmin=0.2 ; b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5050("bruce/rundir5050", "Bruce 5050", "Bruce Shaw", cal(2020, 10, 5),
				"Default5046.   bigger   b= .010",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5051("bruce/rundir5051", "Bruce 5051", "Bruce Shaw", cal(2020, 10, 5),
				"Default5046.   smaller  b= .008",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5052("bruce/rundir5052", "Bruce 5052", "Bruce Shaw", cal(2020, 10, 5),
				"Default5046.   bigger   s2ddf= .9",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5053("bruce/rundir5053", "Bruce 5053", "Bruce Shaw", cal(2020, 10, 5),
				"Default5046.   smaller  s2ddf= .7",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5054("bruce/rundir5054", "Bruce 5054", "Bruce Shaw", cal(2020, 10, 5),
				"Default5046.   bigger   ddfmin= .27",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5055("bruce/rundir5055", "Bruce 5055", "Bruce Shaw", cal(2020, 10, 5),
				"Default5046.   smaller  ddfmin= .23",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5056("bruce/rundir5056", "Bruce 5056", "Bruce Shaw", cal(2020, 10, 5),
				"Default5046.   bigger   tCausalF= .70",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5057("bruce/rundir5057", "Bruce 5057", "Bruce Shaw", cal(2020, 10, 5),
				"Default5046.   smaller  tCausalF= .65",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5133("bruce/rundir5133", "Bruce 5133", "Bruce Shaw", cal(2021, 4, 14),
				"Western US.  dx=2.0 Hmax=18.   a=.001 b=.0065",
				null, null),
		BRUCE_5212("bruce/rundir5212", "Bruce 5212", "Bruce Shaw", cal(2021, 4, 27),
				"CA shallower.  H=16km.  b=.009 a=.001",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5250("bruce/rundir5250", "Bruce 5250", "Bruce Shaw", cal(2021, 10, 27),
				"default=same as 4983;  dtau=.67",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5251("bruce/rundir5251", "Bruce 5251", "Bruce Shaw", cal(2021, 10, 27),
				"default but dtau=.65",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5252("bruce/rundir5252", "Bruce 5252", "Bruce Shaw", cal(2021, 10, 27),
				"default but dtau=.66",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5253("bruce/rundir5253", "Bruce 5253", "Bruce Shaw", cal(2021, 10, 27),
				"default but dtau=.68",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5254("bruce/rundir5254", "Bruce 5254", "Bruce Shaw", cal(2021, 10, 27),
				"default but dtau=.69",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5255("bruce/rundir5255", "Bruce 5255", "Bruce Shaw", cal(2021, 10, 29),
				"default but b=.007",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5256("bruce/rundir5256", "Bruce 5256", "Bruce Shaw", cal(2021, 10, 29),
				"default but b=.008",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5257("bruce/rundir5257", "Bruce 5257", "Bruce Shaw", cal(2021, 10, 29),
				"default but b=.010",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5258("bruce/rundir5258", "Bruce 5258", "Bruce Shaw", cal(2021, 10, 29),
				"default but b=.011",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5260("bruce/rundir5260", "Bruce 5260", "Bruce Shaw", cal(2021, 11, 6),
				"default but b=.008,  dtau=.66",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5261("bruce/rundir5261", "Bruce 5261", "Bruce Shaw", cal(2021, 11, 6),
				"default but b=.008,  dtau=.68",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5262("bruce/rundir5262", "Bruce 5262", "Bruce Shaw", cal(2021, 11, 6),
				"default but b=.010,  dtau=.66",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5263("bruce/rundir5263", "Bruce 5263", "Bruce Shaw", cal(2021, 11, 6),
				"default but b=.010,  dtau=.68",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5285("bruce/rundir5285", "Bruce 5285", "Bruce Shaw", cal(2022, 1, 17),
				"Lower b in shallow:  b_shallow= .003, b_deep=.009, h_shallow=3.0km, otherwise default r4983",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5286("bruce/rundir5286", "Bruce 5286", "Bruce Shaw", cal(2022, 1, 21),
				"Lower b in shallow:  b_shallow= .003,  b_deep=.010, h_shallow=3.0km, otherwise default r4983",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5287("bruce/rundir5287", "Bruce 5287", "Bruce Shaw", cal(2022, 1, 21),
				"Lower b in shallow:  b_shallow= .003,  b_deep=.011, h_shallow=3.0km, otherwise default r4983",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5288("bruce/rundir5288", "Bruce 5288", "Bruce Shaw", cal(2022, 1, 21),
				"Lower b in shallow:  b_shallow= .0025, b_deep=.011, h_shallow=3.0km, otherwise default r4983",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5289("bruce/rundir5289", "Bruce 5289", "Bruce Shaw", cal(2022, 1, 21),
				"Lower b in shallow:  b_shallow= .002,  b_deep=.011,   state2f=.15,  otherwise default r4983",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5292("bruce/rundir5292", "Bruce 5292", "Bruce Shaw", cal(2022, 1, 20),
				"Higher res; b=.010",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5293("bruce/rundir5293", "Bruce 5293", "Bruce Shaw", cal(2022, 1, 21),
				"Higher res:  b_shallow= .003, b_deep=.011, h_shallow=3.0km, otherwise default r4983",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5310("bruce/rundir5310", "Bruce 5310", "Bruce Shaw", cal(2022, 1, 17),
				"Vmax=3.0 m/s; otherwise default r4983",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5321("bruce/rundir5321", "Bruce 5321", "Bruce Shaw", cal(2022, 1, 24),
				"Lower b in shallow: b_shallow= .002, b_deep=.011, state2minf= .05, state2f=.9,",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5294("bruce/rundir5294", "Bruce 5294", "Bruce Shaw", cal(2022, 1, 25),
				"Higher res: b_shallow= .002, b_deep=.013, h_shallow=3.0km, state2minf= .15, state2f=.8",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5320("bruce/rundir5320", "Bruce 5320", "Bruce Shaw", cal(2022, 1, 25),
				"b_shallow= .002, b_deep=.011, h_shallow=3.0km, state2minf= .10, state2f=.8",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5295("bruce/rundir5295", "Bruce 5295", "Bruce Shaw", cal(2022, 1, 30),
				"Higher res: b_shallow= .003, b_deep=.013, h_shallow=3.0km, state2minf= .20, state2f=.8",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5296("bruce/rundir5296", "Bruce 5296", "Bruce Shaw", cal(2022, 1, 30),
				"Higher res: b_shallow= .003, b_deep=.015, h_shallow=3.0km, state2minf= .20, state2f=.8",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5410("bruce/rundir5410", "Bruce 5410", "Bruce Shaw", cal(2022, 3, 4),
				"b_shallow= .003, b_deep=.012, h_shallow=3.0km, state2minf= .3, state2f=.8",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5411("bruce/rundir5411", "Bruce 5411", "Bruce Shaw", cal(2022, 3, 4),
				"b_shallow= .0025, b_deep=.013, h_shallow=3.0km, state2minf= .3, state2f=.8",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5412("bruce/rundir5412", "Bruce 5412", "Bruce Shaw", cal(2022, 3, 5),
				"b_shallow= .0025, b_deep=.012, h_shallow=3.0km, state2minf= .3, state2f=.8",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5413("bruce/rundir5413", "Bruce 5413", "Bruce Shaw", cal(2022, 3, 11),
				"b_shallow= .0025, b_deep=.012, h_shallow=3.0km, state2minf= .2, state2f=.8",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5414("bruce/rundir5414", "Bruce 5414", "Bruce Shaw", cal(2022, 3, 11),
				"b_shallow= .0025, b_deep=.012, h_shallow=3.0km, state2minf= .15, state2f=.8",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5415("bruce/rundir5415", "Bruce 5415", "Bruce Shaw", cal(2022, 3, 11),
				"b_shallow= .0025, b_deep=.012, h_shallow=3.0km, state2minf= .1, state2f=.8",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5416("bruce/rundir5416", "Bruce 5416", "Bruce Shaw", cal(2022, 4, 20),
				"b_shallow= .002, b_deep=.013, h_shallow=3.0km, state2minf= .20, state2f=.8",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5418("bruce/rundir5418", "Bruce 5418", "Bruce Shaw", cal(2022, 4, 20),
				"b_shallow= .002, b_deep=.012, h_shallow=3.0km, state2minf= .20, state2f=.8, tdelay= 0.65",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5442("bruce/rundir5442", "Bruce 5442", "Bruce Shaw", cal(2022, 5, 12),
				"high res delta=1.0km, b_shallow= .002, b_deep=.012, h_shallow=3.0km, state2minf= .20, state2f=.8, tdelay= 0.67",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5443("bruce/rundir5443", "Bruce 5443", "Bruce Shaw", cal(2022, 5, 12),
				"high res delta=1.0km, b_shallow= .005, b_deep=.010, h_shallow=3.0km, state2minf= .20, state2f=.8, tdelay= 0.67",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5444("bruce/rundir5444", "Bruce 5444", "Bruce Shaw", cal(2022, 5, 12),
				"high res delta=1.0km, b_shallow= .002, b_deep=.014, h_shallow=3.0km, state2minf= .20, state2f=.8, tdelay= 0.67",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5445("bruce/rundir5445", "Bruce 5445", "Bruce Shaw", cal(2022, 5, 17),
				"high res delta=1.0km, b_shallow= .003, b_deep=.013, h_shallow=3.0km, state2minf= .20, state2f=.8, tdelay= 0.67",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5446("bruce/rundir5446", "Bruce 5446", "Bruce Shaw", cal(2022, 5, 17),
				"high res delta=1.0km, b_shallow= .003, b_deep=.013, h_shallow=3.0km, state2minf= .20, state2f=.8, tdelay= 0.66",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5447("bruce/rundir5447", "Bruce 5447", "Bruce Shaw", cal(2022, 5, 17),
				"high res delta=1.0km, b_shallow= .003, b_deep=.013, h_shallow=3.0km, state2minf= .20, state2f=.8, tdelay= 0.68",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5448("bruce/rundir5448", "Bruce 5448", "Bruce Shaw", cal(2022, 5, 22),
				"high res delta=1.0km, b_shallow= .003, b_deep=.013, h_shallow=3.0km, state2minf= .20, state2f=.8, tdelay= 0.64",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5449("bruce/rundir5449", "Bruce 5449", "Bruce Shaw", cal(2022, 5, 22),
				"high res delta=1.0km, b_shallow= .003, b_deep=.013, h_shallow=3.0km, state2minf= .20, state2f=.8, tdelay= 0.62",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5450("bruce/rundir5450", "Bruce 5450", "Bruce Shaw", cal(2022, 5, 22),
				"high res delta=1.0km, b_shallow= .004, b_deep=.013, h_shallow=3.0km, state2minf= .20, state2f=.8, tdelay= 0.67",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5552("bruce/rundir5552", "Bruce 5552", "Bruce Shaw", cal(2023, 2, 25),
				"NZ model.  default single layer dynamic friction",
				59, 'G'),
		BRUCE_5554("bruce/rundir5554", "Bruce 5554", "Bruce Shaw", cal(2023, 2, 25),
				"NZ two layer.  b_shallow= .004, b_deep=.013, h_shallow=3.0km, state2minf= .30, state2f=.8, tdelay= 0.67",
				59, 'G'),
		BRUCE_5566("rundir5566", "Bruce 5566", "Bruce Shaw", cal(2023, 3, 25),
				"NZ dynamic,  bdeep=.01, bshallow=.003,  alpha=0.25",
				59, 'G'),
		BRUCE_5566_SUB("rundir5566_subduction", "Bruce 5566 (Subduction Only)", "Bruce Shaw", cal(2023, 3, 25),
				"NZ dynamic,  bdeep=.01, bshallow=.003,  alpha=0.25; Filtered for slip on subduction patches only;",
				59, 'G'),
		BRUCE_5566_CRUSTAL("rundir5566_crustal", "Bruce 5566 (Crustal Only)", "Bruce Shaw", cal(2023, 3, 25),
				"NZ dynamic,  bdeep=.01, bshallow=.003,  alpha=0.25; Filtered for slip on crustal patches only;",
				59, 'G'),
		BRUCE_5585("rundir5585", "Bruce 5585", "Bruce Shaw", cal(2023, 8, 2),
				"NZ dynamic, bdeep=.009, bRampto=.003, alpha=0.25",
				59, 'G'),
		BRUCE_5595("rundir5595", "Bruce 5595", "Bruce Shaw", cal(2023, 8, 17),
				"NZ dynamic, bdeep=.009, bRampto=.002, alpha=0.25, hload=hst=3.0",
				59, 'G'),
		BRUCE_5597("rundir5597", "Bruce 5597", "Bruce Shaw", cal(2023, 8, 17),
				"NZ dynamic, bdeep=.009, bshallow=.003, alpha=0.25, hload=hst=3.0",
				59, 'G'),
		BRUCE_5597_SUB("rundir5597_subduction", "Bruce 5597 (Subduction Only)", "Bruce Shaw", cal(2023, 8, 17),
				"NZ dynamic, bdeep=.009, bshallow=.003, alpha=0.25, hload=hst=3.0; Filtered for slip on subduction patches only;",
				59, 'G'),
		BRUCE_5597_CRUSTAL("rundir5597_crustal", "Bruce 5597 (Crustal Only)", "Bruce Shaw", cal(2023, 8, 17),
				"NZ dynamic, bdeep=.009, bshallow=.003, alpha=0.25, hload=hst=3.0; Filtered for slip on crustal patches only;",
				59, 'G'),
		BRUCE_5615("rundir5615", "Bruce 5615", "Bruce Shaw", cal(2023, 8, 24),
				"CA high res; dynamic, bdeep=.013, bshallow=.003, alpha=0.25, hload=hst=3.0",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5652("rundir5652", "Bruce 5652", "Bruce Shaw", cal(2023, 8, 28),
				"CA high res; dynamic, bdeep=.013, bshallow=.003, alpha=0.25, hload=hst=3.0",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5653("rundir5653", "Bruce 5653", "Bruce Shaw", cal(2023, 8, 28),
				"CA high res; dynamic, bdeep=.012, bshallow=.0025, alpha=0.25, hload=hst=3.0",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5654("rundir5654", "Bruce 5654", "Bruce Shaw", cal(2023, 8, 28),
				"CA high res; dynamic, bdeep=.011, bshallow=.003, alpha=0.25, hload=hst=3.0",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5658("rundir5658", "Bruce 5658", "Bruce Shaw", cal(2023, 9, 20),
				"CA high res; dynamic, bdeep=.012, bshallow=.003, alpha=0.25, hload=hst=3.0, fractionLow=.10, blow=.003, Uniform backslip non-hybrid loaded",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5671("rundir5671", "Bruce 5671", "Bruce Shaw", cal(2023, 11, 6),
				"CA high res; dynamic, bdeep=.012, bshallow=.003, alpha=0.25, hload=hst=3.0, fractionLow=.10, blow=.003",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5672("rundir5672", "Bruce 5672", "Bruce Shaw", cal(2023, 11, 6),
				"CA high res; dynamic, bdeep=.014, bshallow=.003, alpha=0.25, hload=hst=3.0, fractionLow=.10, blow=.003",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_5684("rundir5684", "Bruce 5684", "Bruce Shaw", cal(2023, 12, 7),
				"NZ dynamic, bdeep=.010, bshallow=.003, alpha=0.25, hload=hst=3.0",
				59, 'G'),
		BRUCE_5684_SUB("rundir5684_subduction", "Bruce 5684 (Subduction Only)", "Bruce Shaw", cal(2023, 12, 7),
				"NZ dynamic, bdeep=.010, bshallow=.003, alpha=0.25, hload=hst=3.0; Filtered for slip on subduction patches only;",
				59, 'G'),
		BRUCE_5684_CRUSTAL("rundir5684_crustal", "Bruce 5684 (Crustal Only)", "Bruce Shaw", cal(2023, 12, 7),
				"NZ dynamic, bdeep=.010, bshallow=.003, alpha=0.25, hload=hst=3.0; Filtered for slip on crustal patches only;",
				59, 'G'),
		BRUCE_5685("rundir5685", "Bruce 5685", "Bruce Shaw", cal(2023, 12, 7),
				"NZ dynamic, bdeep=.010, bshallow=.003, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.015",
				59, 'G'),
		BRUCE_5685_SUB("rundir5685_subduction", "Bruce 5685 (Subduction Only)", "Bruce Shaw", cal(2023, 12, 7),
				"NZ dynamic, bdeep=.010, bshallow=.003, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.015; Filtered for slip on subduction patches only;",
				59, 'G'),
		BRUCE_5685_CRUSTAL("rundir5685_crustal", "Bruce 5685 (Crustal Only)", "Bruce Shaw", cal(2023, 12, 7),
				"NZ dynamic, bdeep=.010, bshallow=.003, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.015; Filtered for slip on crustal patches only;",
				59, 'G'),
		BRUCE_5687("rundir5687", "Bruce 5687", "Bruce Shaw", cal(2023, 12, 26),
				"NZ dynamic, bdeep=.009, bshallow=.003, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.018",
				59, 'G'),
		BRUCE_5687_SUB("rundir5687_subduction", "Bruce 5687 (Subduction Only)", "Bruce Shaw", cal(2023, 12, 26),
				"NZ dynamic, bdeep=.009, bshallow=.003, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.018; Filtered for slip on subduction patches only;",
				59, 'G'),
		BRUCE_5687_CRUSTAL("rundir5687_crustal", "Bruce 5687 (Crustal Only)", "Bruce Shaw", cal(2023, 12, 26),
				"NZ dynamic, bdeep=.009, bshallow=.003, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.018; Filtered for slip on crustal patches only;",
				59, 'G'),
		BRUCE_5689("rundir5689", "Bruce 5689", "Bruce Shaw", cal(2023, 12, 26),
				"NZ dynamic, bdeep=.009, bshallow=.002, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.018",
				59, 'G'),
		BRUCE_5689_SUB("rundir5689_subduction", "Bruce 5689 (Subduction Only)", "Bruce Shaw", cal(2023, 12, 26),
				"NZ dynamic, bdeep=.009, bshallow=.002, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.018; Filtered for slip on subduction patches only;",
				59, 'G'),
		BRUCE_5689_CRUSTAL("rundir5689_crustal", "Bruce 5689 (Crustal Only)", "Bruce Shaw", cal(2023, 12, 26),
				"NZ dynamic, bdeep=.009, bshallow=.002, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.018; Filtered for slip on crustal patches only;",
				59, 'G'),
		BRUCE_5691("rundir5691", "Bruce 5691", "Bruce Shaw", cal(2024, 1, 5),
				"NZ dynamic, bdeep=.009, bshallow=.002, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.019",
				59, 'G'),
		BRUCE_5691_SUB("rundir5691_subduction", "Bruce 5691 (Subduction Only)", "Bruce Shaw", cal(2024, 1, 5),
				"5691, but filtered for slip on subduction patches only;",
				59, 'G'),
		BRUCE_5691_CRUSTAL("rundir5691_crustal", "Bruce 5691 (Crustal Only)", "Bruce Shaw", cal(2024, 1, 5),
				"5691, but filtered for slip on crustal patches only;",
				59, 'G'),
		BRUCE_5696("rundir5696", "Bruce 5696", "Bruce Shaw", cal(2024, 1, 19),
				"NZ dynamic, bdeep=.010, bshallow=.002, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.018",
				59, 'G'),
		BRUCE_5696_SUB("rundir5696_subduction", "Bruce 5696 (Subduction Only)", "Bruce Shaw", cal(2024, 1, 19),
				"5696, but filtered for slip on subduction patches only;",
				59, 'G'),
		BRUCE_5696_CRUSTAL("rundir5696_crustal", "Bruce 5696 (Crustal Only)", "Bruce Shaw", cal(2024, 1, 19),
				"5696, but filtered for slip on crustal patches only;",
				59, 'G'),
		BRUCE_5697("rundir5697", "Bruce 5697", "Bruce Shaw", cal(2024, 1, 19),
				"NZ dynamic, bdeep=.010, bshallow=.0015, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.018",
				59, 'G'),
		BRUCE_5697_SUB("rundir5697_subduction", "Bruce 5697 (Subduction Only)", "Bruce Shaw", cal(2024, 1, 19),
				"5697, but filtered for slip on subduction patches only;",
				59, 'G'),
		BRUCE_5697_CRUSTAL("rundir5697_crustal", "Bruce 5697 (Crustal Only)", "Bruce Shaw", cal(2024, 1, 19),
				"5697, but filtered for slip on crustal patches only;",
				59, 'G'),
		BRUCE_5774("rundir5774", "Bruce 5774", "Bruce Shaw", cal(2024, 3, 5),
				"NZ dynamic, sigma0=200, bdeep=.005, bshallow=.003, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.012",
				59, 'G'),
		BRUCE_5774_SUB("rundir5774_subduction", "Bruce 5774 (Subduction Only)", "Bruce Shaw", cal(2024, 3, 5),
				"5774, but filtered for slip on subduction patches only;",
				59, 'G'),
		BRUCE_5774_CRUSTAL("rundir5774_crustal", "Bruce 5774 (Crustal Only)", "Bruce Shaw", cal(2024, 3, 5),
				"5774, but filtered for slip on crustal patches only;",
				59, 'G'),
		BRUCE_5775("rundir5775", "Bruce 5775", "Bruce Shaw", cal(2024, 3, 5),
				"NZ dynamic, sigma0=200, bdeep=.005, bshallow=.002, alpha=0.25, hload=hst=3.0, fractionLow=.10, bhigh=.012",
				59, 'G'),
		BRUCE_5775_SUB("rundir5775_subduction", "Bruce 5775 (Subduction Only)", "Bruce Shaw", cal(2024, 3, 5),
				"5775, but filtered for slip on subduction patches only;",
				59, 'G'),
		BRUCE_5775_CRUSTAL("rundir5775_crustal", "Bruce 5775 (Crustal Only)", "Bruce Shaw", cal(2024, 3, 5),
				"5775, but filtered for slip on crustal patches only;",
				59, 'G'),
		BRUCE_5844("rundir5892", "Bruce 5892", "Bruce Shaw", cal(2024, 8, 25),
				"WUSav, delta=2.0km, sigma0=100, b=.008, alpha=0.25",
				NSHM23_FaultModels.WUS_FM_v3, NSHM23_DeformationModels.AVERAGE);
		
		private String dirName;
		private RSQSimCatalog catalog;
		
		private Catalogs(String dirName, String name, String author, GregorianCalendar date, String metadata,
				RupSetFaultModel fm, RupSetDeformationModel dm) {
			this(dirName, name, author, date, metadata, fm, dm, -1, '0');
		}
		
		private Catalogs(String dirName, String name, String author, GregorianCalendar date, String metadata,
				int utmZone, char utmBand) {
			this(dirName, name, author, date, metadata, null, null, utmZone, utmBand);
		}
		
		private Catalogs(String dirName, String name, String author, GregorianCalendar date, String metadata,
				RupSetFaultModel fm, RupSetDeformationModel dm, int utmZone, char utmBand) {
			this.dirName = dirName;
			if (utmZone < 0) {
				if (fm != null && fm instanceof NZ_CompFaultModel) {
					utmZone = 59;
					utmBand = 'G';
				} else {
					utmZone = 11;
					utmBand = 'N';
				}
			}
			catalog = new RSQSimCatalog(name, author, date, metadata, fm, dm, utmZone, utmBand);
		}
		
		public RSQSimCatalog instance() {
			for (File baseDir : catalogLocations) {
				File dir = new File(baseDir, dirName);
				if (dir.exists()) {
					catalog.dir = dir;
					return catalog;
				}
			}
			throw new IllegalStateException("Couldn't locate catalog "+name());
		}
		
		public RSQSimCatalog instance(File baseDir) {
			File dir = new File(baseDir, dirName);
			Preconditions.checkState(dir.exists(), "Catalog dir doesn't exist: %s", dir.getAbsolutePath());
			catalog.dir = dir;
			return catalog;
		}
	}
	
	private File dir;
	private String name;
	private String author;
	private GregorianCalendar date;
	private String metadata;
	private RupSetFaultModel fm;
	private RupSetDeformationModel dm;
	
	private int utmZone;
	private char utmBand;
	
	private double constSlipVel = Double.NaN;
	private Map<Integer, Double> slipVels = null;
	private double aveArea = Double.NaN;
	private int numEvents = -1;
	private double durationYears = Double.NaN;
	private Map<String, String> params;
	private Boolean variableSlipSpeed = null;
	
	private List<SimulatorElement> elements;
	private RSQSimStateTransitionFileReader transReader;
	private List<? extends FaultSection> subSects;
	private Map<Integer, Double> subSectAreas;
	private Map<IDPairing, Double> subSectDistsCache;
	private RSQSimSubSectionMapper subSectMapper;
	
	private static final File fmDmSolDir = new File(System.getProperty("user.home"), ".opensha/ucerf3_fm_dm_sols/");
	private FaultSystemSolution compSol;
	
	public static final double MIN_SUB_SECT_FRACT_DEFAULT = 0.2;
	
	private double minFractForInclusion = MIN_SUB_SECT_FRACT_DEFAULT;
	
	private static Table<RupSetFaultModel, RupSetDeformationModel, FaultSystemSolution> compSolsTable = HashBasedTable.create();
	
	public static final String XML_METADATA_NAME = "RSQSimCatalog";
	
	private RSQSimCatalog(String name, String author, GregorianCalendar date, String metadata,
			RupSetFaultModel fm, RupSetDeformationModel dm, int utmZone, char utmBand) {
		this(null, name, author, date, metadata, fm, dm, utmZone, utmBand);
	}
	
	/**
	 * Creates a new RSQSimCatalog instance without any the metadata
	 * @param dir directory containing the list/geometry files
	 * @param name name of this catalog
	 * @param fm fault model for comparisons
	 * @param dm defomration model for comparisons
	 */
	public RSQSimCatalog(File dir, String name, RupSetFaultModel fm, RupSetDeformationModel dm) {
		this(dir, name, null, null, null, fm, dm);
	}

	/**
	 * Creates a new RSQSimCatalog instance with all metadata
	 * 
	 * @param dir directory containing the list/geometry files
	 * @param name name of this catalog
	 * @param author author of this catalog
	 * @param date creation date for this catalog
	 * @param metadata description of this catalog
	 * @param fm fault model for comparisons
	 * @param dm defomration model for comparisons
	 */
	public RSQSimCatalog(File dir, String name, String author, GregorianCalendar date, String metadata,
			RupSetFaultModel fm, RupSetDeformationModel dm) {
		this(dir, name, author, date, metadata, fm, dm, 11, 'N');
	}

	/**
	 * Creates a new RSQSimCatalog instance with all metadata
	 * 
	 * @param dir directory containing the list/geometry files
	 * @param name name of this catalog
	 * @param author author of this catalog
	 * @param date creation date for this catalog
	 * @param metadata description of this catalog
	 * @param fm fault model for comparisons
	 * @param dm defomration model for comparisons
	 */
	public RSQSimCatalog(File dir, String name, String author, GregorianCalendar date, String metadata,
			RupSetFaultModel fm, RupSetDeformationModel dm, int utmZone, char utmBand) {
		this.dir = dir;
		this.name = name;
		this.author = author;
		this.date = date;
		this.metadata = metadata;
		this.fm = fm;
		this.dm = dm;
		this.utmZone = utmZone;
		this.utmBand = utmBand;
	}
	
	public File getCatalogDir() {
		return dir;
	}

	public String getName() {
		return name;
	}

	public String getAuthor() {
		return author;
	}

	public GregorianCalendar getDate() {
		return date;
	}

	public String getMetadata() {
		return metadata;
	}
	
	public RupSetFaultModel getFaultModel() {
		return fm;
	}

	public RupSetDeformationModel getDeformationModel() {
		return dm;
	}
	
	public synchronized boolean isVariableSlipSpeed() throws IOException {
		if (variableSlipSpeed != null)
			return variableSlipSpeed;
		Map<String, String> params = getParams();
		
		String variableSlipSpeedStr = params.get("variableSlipSpeed");
		if (variableSlipSpeedStr  == null) {
			variableSlipSpeed = false;
		} else {
			int variableSlipSpeedInt = Integer.parseInt(variableSlipSpeedStr);
			variableSlipSpeed = variableSlipSpeedInt > 0;
		}
		
		return variableSlipSpeed;
	}
	
	public synchronized Map<Integer, Double> getSlipVelocities() throws IOException {
		if (slipVels == null) {
			if (Doubles.isFinite(constSlipVel) && constSlipVel > 0) {
				List<SimulatorElement> elems = getElements();
				slipVels = new HashMap<>();
				for (int i=0; i<elems.size(); i++)
					slipVels.put(elems.get(i).getID(), constSlipVel);
			} else {
				Map<String, String> params = getParams();
				Preconditions.checkNotNull(params, "Params not found");
				String ddotEQFname = params.get("ddotEQFname");
				if (ddotEQFname != null && !ddotEQFname.trim().isEmpty()) {
					File ddotEQFile = new File(getCatalogDir(), ddotEQFname);
					Preconditions.checkState(ddotEQFile.exists(),
							"ddotEQFname = %s doesn't exist in %s", ddotEQFname, getCatalogDir().getAbsolutePath());
					double[] velArray = loadDoubleInputFile(ddotEQFile);
					List<SimulatorElement> elems = getElements();
					Preconditions.checkState(velArray.length == elems.size(), "expected %s patch velocities, have %s", elems.size(), velArray.length);
					slipVels = new HashMap<>();
					for (int i=0; i<elems.size(); i++)
						slipVels.put(elems.get(i).getID(), velArray[i]);
					constSlipVel = Double.NaN;
				} else {
					String ddotEQ = params.get("ddotEQ_1");
					Preconditions.checkNotNull(ddotEQ, "ddotEQ_1 not in params file");
					constSlipVel = Double.parseDouble(ddotEQ);
					List<SimulatorElement> elems = getElements();
					slipVels = new HashMap<>();
					for (int i=0; i<elems.size(); i++)
						slipVels.put(elems.get(i).getID(), constSlipVel);
				}
			}
		}
		return slipVels;
	}
	
	private static double[] loadDoubleInputFile(File file) throws IOException {
		List<Double> vals = new ArrayList<>();
		for (String line : Files.readLines(file, Charset.defaultCharset())) {
			line = line.trim();
			if (line.isEmpty() || line.startsWith("#"))
				continue;
			vals.add(Double.parseDouble(line));
		}
		return Doubles.toArray(vals);
	}
	
	private static DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd");
	private static DecimalFormat areaDF = new DecimalFormat("0.00");
	private static DecimalFormat groupedIntDF = new DecimalFormat("#");
	static {
		groupedIntDF.setGroupingUsed(true);
		groupedIntDF.setGroupingSize(3);
	}
	
	public List<String> getMarkdownMetadataTable() {
		TableBuilder builder = MarkdownUtils.tableBuilder();
		builder.addLine("**Catalog**", getName());
		String authorLine = getAuthor() == null ? "(unknown)" : getAuthor();
		if (getDate() != null)
			authorLine += ", "+dateFormat.format(getDate().getTime());
		builder.addLine("**Author**", authorLine);
		builder.addLine("**Description**", getMetadata());
		builder.addLine("**Fault/Def Model**", fm+", "+dm);
		try {
			Map<Integer, Double> slipVels = getSlipVelocities();
			String velStr;
			if (isVariableSlipSpeed()) {
				velStr = "Fully Variable";
			} else if (Double.isFinite(constSlipVel)) {
				velStr = (float)constSlipVel+" m/s";
			} else {
				MinMaxAveTracker velTrack = new MinMaxAveTracker();
				for (Double vel : slipVels.values())
					velTrack.addValue(vel);
				velStr = "Patch Variable, range=["+(float)velTrack.getMin()+" "+(float)velTrack.getMax()+"], mean="+(float)velTrack.getAverage();
			}
			builder.addLine("**Slip Velocity**", velStr);
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			double  aveArea = getAveArea();
			builder.addLine("**Average Element Area**", areaDF.format(aveArea)+" km^2");
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			int numEvents = getNumEvents();
			double durationYears = getDurationYears();
			builder.addLine("**Length**", groupedIntDF.format(numEvents)+" events in "
					+groupedIntDF.format(durationYears)+" years");
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			Map<String, String> params = getParams();
			if (params != null) {
				double a = Double.parseDouble(params.get("A_1"));
				double b = Double.parseDouble(params.get("B_1"));
				String ddotEQ = params.get("ddotEQ_1");
				if (ddotEQ.contains("."))
					ddotEQ = Float.parseFloat(ddotEQ)+"";
				builder.addLine("**Frictional Params**", "a="+(float)a+", b="+(float)b+", (b-a)="+(float)(b-a)+", ddotEQ="+ddotEQ);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return builder.build();
	}
	
	public synchronized double getAveArea() throws IOException {
		if (Double.isNaN(aveArea)) {
			List<SimulatorElement> elements = getElements();
			aveArea = 0d;
			for (SimulatorElement e : elements)
				aveArea += e.getArea()*1e-6;
			aveArea /= elements.size();
		}
		return aveArea;
	}
	
	public synchronized int getNumEvents() throws IOException {
		if (numEvents < 0)
			numEvents = RSQSimFileReader.getNumEvents(getCatalogDir());
		return numEvents;
	}
		
	public synchronized double getDurationYears() throws IOException {
		if (Double.isNaN(durationYears))
			durationYears = RSQSimFileReader.getDurationYears(getCatalogDir());
		return durationYears;
	}
	
	public synchronized Map<String, String> getParams() throws IOException {
		if (params == null) {
			File paramFile = getParamFile();
			if (paramFile != null)
				params = RSQSimFileReader.readParams(paramFile);
		}
		return params;
	}
	
	public File getParamFile() throws IOException {
		return RSQSimFileReader.getParamFile(getCatalogDir());
	}
	
	public void writeMarkdownSummary(File dir, boolean plots, boolean replot,
			StandardPlots... standardPlots) throws IOException {
		writeMarkdownSummary(dir, plots, replot, 6d, standardPlots);
	}
	
	public void writeMarkdownSummary(File dir, boolean plots, boolean replot, double plotMinMag,
			StandardPlots... standardPlots) throws IOException {
		List<String> lines = new LinkedList<>();
		String topLink = "*[(top)](#"+MarkdownUtils.getAnchorName(getName())+")*";
		lines.add("# "+getName());
		lines.add("## Metadata");
		lines.addAll(getMarkdownMetadataTable());
		lines.add("");
		int tocIndex = lines.size();
		
		List<String> hazardLinks = new ArrayList<>();
		List<String> hazardNames = new ArrayList<>();
		
		List<String> hazardClusterLinks = new ArrayList<>();
		List<String> hazardClusterNames = new ArrayList<>();
		
		List<String> occCopulaLinks = new ArrayList<>();
		List<String> occCopulaNames = new ArrayList<>();
		
		List<String> multiFaultLinks = new ArrayList<>();
		List<String> multiFaultNames = new ArrayList<>();

		String partBSummaryLink = null;
		String crustalSubductionGMsLink = null;
		String vmCompareRotRupLink = null;
		String extremeEventLink = null;
		String parentMFDLink = null;
		String distCalcLink = null;
		String searchEventLink = null;
		
		File[] dirList = dir.listFiles();
		Arrays.sort(dirList, new FileNameComparator());
		for (File subDir : dirList) {
			if (!subDir.isDirectory())
				continue;
			File mdFile = new File(subDir, "README.md");
			if (!mdFile.exists())
				continue;
			String name = subDir.getName();
			if (name.startsWith("hazard_")) {
				String hazName;
				if (name.contains("_pga")) {
					hazName = "PGA";
				} else if (name.contains("t_dependence")) {
					hazName = "T Dependence";
				} else {
					Preconditions.checkState(name.contains("_sa"));
					String periodStr = name.substring(name.indexOf("_sa_")+4);
					periodStr = periodStr.substring(0, periodStr.indexOf("s"));
					double period = Double.parseDouble(periodStr);
					hazName = (float)period+"s SA";
				}
				
				if (name.contains("_sigma")) {
					String sigmaStr = name.substring(name.indexOf("_sigma")+6);
					if (sigmaStr.contains("_"))
						sigmaStr = sigmaStr.substring(0, sigmaStr.indexOf("_"));
					hazName += ", Fixed σ="+sigmaStr;
				}
				
				if (name.contains("_gmpe")) {
					String gmpeStr = name.substring(name.indexOf("_gmpe")+5);
					if (gmpeStr.contains("_"))
						gmpeStr = gmpeStr.substring(0, gmpeStr.indexOf("_"));
					hazName += ", "+gmpeStr;
				}
				
				if (name.contains("_sectArea")) {
					String areaStr = name.substring(name.indexOf("_sectArea")+9);
					if (areaStr.contains("_"))
						areaStr = areaStr.substring(0, areaStr.indexOf("_"));
					hazName += ", SectAreaFract="+areaStr;
				}
				
				if (name.contains("_cluster")) {
					hazardClusterLinks.add(name);
					hazardClusterNames.add(hazName);
				} else {
					hazardLinks.add(name);
					hazardNames.add(hazName);
				}
			} else if (name.startsWith("multi_fault")) {
				if (name.equals("multi_fault")) {
					multiFaultLinks.add(name);
					multiFaultNames.add("Default sub-section mapping");
				} else {
					Preconditions.checkState(name.contains("_area_fract_"));
					multiFaultLinks.add(name);
					double fract = Double.parseDouble(name.substring(
							name.indexOf("_area_fract_")+("_area_fract_").length()));
					multiFaultNames.add("Sub-section mapping with areaFract="+(fract));
				}
			} else if (name.equals("extreme_events")) {
				Preconditions.checkState(extremeEventLink == null,
						"Duplicate Extreme Event dirs! %s and %s", name, extremeEventLink);
				extremeEventLink = name;
			} else if (name.startsWith("occupancy_copula_m")) {
				String title = MarkdownUtils.getTitle(mdFile);
				occCopulaLinks.add(name);
				occCopulaNames.add(title);
			} else if (name.equals("bbp_vm_rot_rup_compare")) {
				vmCompareRotRupLink = name;
			} else if (name.equals("bbp_part_b_summary")) {
				partBSummaryLink = name;
			} else if (name.equals("crustal_subduction_gms")) {
				crustalSubductionGMsLink = name;
			} else if (name.equals("parent_sect_mfds")) {
				parentMFDLink = name;
			} else if (name.equals("dist_method_comparisons")) {
				distCalcLink = name;
			} else if (name.equals("search_events")) {
				searchEventLink = name;
			}
		}
		
		if (!hazardLinks.isEmpty()) {
			lines.add("");
			lines.add("## Hazard Comparisons");
			lines.add(topLink);
			lines.add("");
			for (int i=0; i<hazardNames.size(); i++)
				lines.add("* ["+hazardNames.get(i)+"]("+hazardLinks.get(i)+"/)");
		}
		if (!hazardClusterLinks.isEmpty()) {
			lines.add("");
			lines.add("## Hazard Clustering Comparisons");
			lines.add(topLink);
			lines.add("");
			for (int i=0; i<hazardClusterNames.size(); i++)
				lines.add("* ["+hazardClusterNames.get(i)+"]("+hazardClusterLinks.get(i)+"/)");
		}
		if (!multiFaultLinks.isEmpty()) {
			lines.add("");
			lines.add("## Multi-Fault Rupture Comparisons");
			lines.add(topLink);
			lines.add("");
			if (multiFaultLinks.size() == 1) {
				lines.add("[Multi-Fault Rupture Comparisons here]("+multiFaultLinks.get(0)+"/)");
			} else {
				for (int i=0; i<multiFaultLinks.size(); i++)
					lines.add("* ["+multiFaultNames.get(i)+"]("+multiFaultLinks.get(i)+"/)");
			}
		}
		if (parentMFDLink != null) {
			lines.add("");
			lines.add("## Parent Section MFDs");
			lines.add(topLink);
			lines.add("");
			lines.add("[Parent Section MFDs here]("+parentMFDLink+"/)");
		}
		if (searchEventLink != null) {
			lines.add("");
			lines.add("## Event Search Results");
			lines.add(topLink);
			lines.add("");
			lines.add("[Event search results here]("+searchEventLink+"/)");
		}
		if (distCalcLink != null) {
			lines.add("");
			lines.add("## Distance Calculation Method Comparisons");
			lines.add(topLink);
			lines.add("");
			lines.add("[Distance Calculation Methods here]("+distCalcLink+"/)");
		}
		if (extremeEventLink != null) {
			lines.add("");
			lines.add("## Extreme Event Examples");
			lines.add(topLink);
			lines.add("");
			lines.add("[Extreme Event Examples Here]("+extremeEventLink+"/)");
		}
		if (!occCopulaLinks.isEmpty()) {
			lines.add("");
			lines.add("## Occupancy Copulas");
			lines.add(topLink);
			lines.add("");
			for (int i=0; i<occCopulaLinks.size(); i++)
				lines.add("* ["+occCopulaNames.get(i)+"]("+occCopulaLinks.get(i)+"/)");
		}
		
		// BBP Pages
		for (VelocityModel vm : VelocityModel.values()) {
			File bbpDir = new File(dir, "bbp_"+vm.name());
			if (!bbpDir.exists())
				continue;
			
			lines.add("");
			lines.add("## BBP Calculations, "+vm+" Velocity Model");
			lines.add(topLink); lines.add("");
			
			List<String> eventLinks = new ArrayList<>();
			List<String> eventNames = new ArrayList<>();
			
			List<String> gmpeLinks = new ArrayList<>();
			List<String> gmpeNames = new ArrayList<>();
			
			List<String> gmpeGriddedLinks = new ArrayList<>();
			List<String> gmpeGriddedNames = new ArrayList<>();
			
			List<String> gmpeRGLinks = new ArrayList<>();
			List<String> gmpeRGNames = new ArrayList<>();
			
			List<String> gmpeNonErgodicLinks = new ArrayList<>();
			List<String> gmpeNonErgodicNames = new ArrayList<>();
			
			List<String> rotatedRupLinks = new ArrayList<>();
			List<String> rotatedRupNames = new ArrayList<>();
			
			List<String> azimuthalRupLinks = new ArrayList<>();
			List<String> azimuthalRupNames = new ArrayList<>();
			
			List<String> partBLinks = new ArrayList<>();
			List<String> partBNames = new ArrayList<>();
			
			List<String> sourceDetailLinks = new ArrayList<>();
			List<String> sourceDetailNames = new ArrayList<>();
			
			Map<String, List<String>> siteHazardLinksTable = new HashMap<>();
			Map<String, List<String>> siteHazardNamesTable = new HashMap<>();
			
			String rotDDLink = null;
			String sourceSiteLink = null;
			
			dirList = bbpDir.listFiles();
			Arrays.sort(dirList, new FileNameComparator());
			for (File subDir : dirList) {
				if (!subDir.isDirectory())
					continue;
				File mdFile = new File(subDir, "README.md");
				if (!mdFile.exists())
					continue;
				String name = subDir.getName();
				if (name.startsWith("event_")) {
					eventNames.add(MarkdownUtils.getTitle(mdFile));
					eventLinks.add(name);
				} else if (name.startsWith("gmpe_bbp_comparisons_")) {
					String subName = name.substring("gmpe_bbp_comparisons_".length());
					
					boolean gridded = name.contains("_GriddedSites");
					if (gridded)
						subName = subName.replaceAll("_GriddedSites", "");
					
					if (subName.contains("_timeScale")) {
						subName = subName.replaceAll("_timeScale", ", Time Scale Factor: ");
						if (subName.contains("_velScale"))
							subName = subName.replaceAll("_velScale", ", Velocities Scaled");
					}
					
					if (subName.contains("_mech_")) {
						subName = subName.replaceAll("_mech_", ", Focal Mechanism: ");
						subName = subName.replaceAll("vert_ss", "Vertical Strike-Slip");
						subName = subName.replaceAll("ss", "Strike-Slip");
						subName = subName.replaceAll("reverse", "Reverse");
						subName = subName.replaceAll("normal", "Normal");
					}
					
					if (gridded) {
						gmpeGriddedNames.add(subName);
						gmpeGriddedLinks.add(name);
					} else {
						gmpeNames.add(subName);
						gmpeLinks.add(name);
					}
				} else if (name.startsWith("gmpe_bbp_non_ergodic_maps_")) {
					String subName = name.substring("gmpe_bbp_non_ergodic_maps_".length());
					gmpeNonErgodicNames.add(subName);
					gmpeNonErgodicLinks.add(name);
				} else if (name.startsWith("gmpe_bbp_rg_comparisons_")) {
					gmpeRGNames.add(name.substring("gmpe_bbp_rg_comparisons_".length()));
					gmpeRGLinks.add(name);
				} else if (name.equals("catalog_rotd_ratio_comparisons")) {
					Preconditions.checkState(rotDDLink == null, "Duplicate RotDD dirs! %s and %s", name, rotDDLink);
					rotDDLink = name;
				} else if (name.equals("source_site_comparisons")) {
					Preconditions.checkState(sourceSiteLink == null, "Duplicate Source/Site dirs! %s and %s", name, sourceSiteLink);
					sourceSiteLink = name;
				} else if (name.startsWith("source_site_")) {
					sourceDetailLinks.add(name);
					sourceDetailNames.add(MarkdownUtils.getTitle(mdFile));
				} else if (name.startsWith("bbp_part_b")) {
					FilterMethod filter = FilterMethod.fromDirName(name);
					name = checkMoveFilter(subDir, filter);
					partBLinks.add(name);
					partBNames.add(filter.getName());
				} else if (name.startsWith("rotated_ruptures_")) {
					FilterMethod filter = FilterMethod.fromDirName(name);
					name = checkMoveFilter(subDir, filter);
					for (Scenario scenario : Scenario.values()) {
						String dirName = "rotated_ruptures_"+scenario.getPrefix();
						if (name.startsWith(dirName+"_timeScale")) {
							String scenName = scenario.getName();
							String timeScale = name.substring(name.indexOf("_timeScale"));
							if (timeScale.contains("_filter"))
								timeScale = timeScale.substring(0, timeScale.indexOf("_filter"));
							timeScale = timeScale.replaceAll("_timeScale", ", Time Scale Factor: ");
							if (timeScale.contains("_velScale"))
								timeScale = timeScale.replaceAll("_velScale", ", Velocities Scaled");
							rotatedRupLinks.add(name);
							rotatedRupNames.add(scenName+timeScale+", "+filter.getName());
						} else if (name.equals(dirName) || name.startsWith(dirName+"_filter")) {
							rotatedRupLinks.add(name);
							rotatedRupNames.add(scenario.getName()+", "+filter.getName());
						}
					}
					for (RuptureType rupType : RuptureType.values()) {
						if (name.equals("rotated_ruptures_mag_dist_"+rupType.getPrefix())) {
							rotatedRupLinks.add(name);
							rotatedRupNames.add(rupType.getName()+", Mag-Dist Bins");
						}
					}
				} else if (name.startsWith("azimuthal_")) {
					for (Scenario scenario : Scenario.values()) {
						String dirName = "azimuthal_"+scenario.getPrefix();
						if (name.startsWith(dirName+"_timeScale")) {
							String scenName = scenario.getName();
							String timeScale = name.substring(name.indexOf("_timeScale"));
							timeScale = timeScale.replaceAll("_timeScale", ", Time Scale Factor: ");
							if (timeScale.contains("_velScale"))
								timeScale = timeScale.replaceAll("_velScale", ", Velocities Scaled");
							azimuthalRupLinks.add(name);
							azimuthalRupNames.add(scenName+timeScale);
						} else if (name.equals(dirName)) {
							azimuthalRupLinks.add(name);
							azimuthalRupNames.add(scenario.getName());
						}
					}
				} else if (name.startsWith("site_hazard_")) {
					String siteName = name.substring("site_hazard_".length());
					siteName = siteName.substring(0, siteName.indexOf("_"));
					String gmpeName = name.substring("site_hazard_".length()+siteName.length()+1);
					
					if (!siteHazardNamesTable.containsKey(gmpeName)) {
						siteHazardLinksTable.put(gmpeName, new ArrayList<>());
						siteHazardNamesTable.put(gmpeName, new ArrayList<>());
					}
					
					siteHazardLinksTable.get(gmpeName).add(name);
					siteHazardNamesTable.get(gmpeName).add(siteName);
				}
			}
			
			if (!eventNames.isEmpty()) {
				lines.add("");
				lines.add("### Single Event Comparisons, "+vm);
				lines.add(topLink);
				lines.add("");
				for (int i=0; i<eventNames.size(); i++)
					lines.add("* ["+eventNames.get(i)+"]("+bbpDir.getName()+"/"+eventLinks.get(i)+"/)");
			}
			if (!gmpeNames.isEmpty() || !gmpeGriddedNames.isEmpty()) {
				lines.add("");
				lines.add("### Full Catalog GMPE Comparisons, "+vm);
				lines.add(topLink);
				lines.add("");
//				System.out.print("Have "+gmpeNames.size()+" regulars and "+gmpeGriddedNames.size()+" gridded");
				boolean both = !gmpeNames.isEmpty() && !gmpeGriddedNames.isEmpty();
				if (both) {
					lines.add("#### Points Of Interest, "+vm);
					lines.add("");
					for (int i=0; i<gmpeNames.size(); i++)
						lines.add("* ["+gmpeNames.get(i)+"]("+bbpDir.getName()+"/"+gmpeLinks.get(i)+"/)");
					lines.add("");
					lines.add("#### Gridded Sites, "+vm);
					lines.add("");
					for (int i=0; i<gmpeGriddedNames.size(); i++)
						lines.add("* ["+gmpeGriddedNames.get(i)+"]("+bbpDir.getName()+"/"+gmpeGriddedLinks.get(i)+"/)");
				} else {
					for (int i=0; i<gmpeNames.size(); i++)
						lines.add("* ["+gmpeNames.get(i)+"]("+bbpDir.getName()+"/"+gmpeLinks.get(i)+"/)");
					for (int i=0; i<gmpeGriddedNames.size(); i++)
						lines.add("* ["+gmpeGriddedNames.get(i)+"]("+bbpDir.getName()+"/"+gmpeGriddedLinks.get(i)+"/)");
				}
			}
			if (!gmpeRGNames.isEmpty()) {
				lines.add("");
				lines.add("### Full Catalog GMPE Comparisons with BBP Rupture Generator, "+vm);
				lines.add(topLink);
				lines.add("");
				for (int i=0; i<gmpeRGNames.size(); i++)
					lines.add("* ["+gmpeRGNames.get(i)+"]("+bbpDir.getName()+"/"+gmpeRGLinks.get(i)+"/)");
			}
			if (!gmpeNonErgodicNames.isEmpty()) {
				lines.add("");
				lines.add("### Full Catalog GMPE Non-Ergodic Map Comparisons, "+vm);
				lines.add(topLink);
				lines.add("");
				for (int i=0; i<gmpeNonErgodicNames.size(); i++)
					lines.add("* ["+gmpeNonErgodicNames.get(i)+"]("+bbpDir.getName()+"/"+gmpeNonErgodicLinks.get(i)+"/)");
			}
			if (rotDDLink != null) {
				lines.add("");
				lines.add("### Full Catalog RotD100/RotD50 Ratios, "+vm);
				lines.add(topLink);
				lines.add("");
				lines.add("[Full Catalog RotD100/RotD50 Ratios Plotted Here]("+bbpDir.getName()+"/"+rotDDLink+"/)");
			}
			if (sourceSiteLink != null || !sourceDetailLinks.isEmpty()) {
				lines.add("");
				lines.add("### Source/Site Ground Motion Comparisons, "+vm);
				lines.add(topLink);
				if (sourceSiteLink != null) {
					lines.add("");
					lines.add("[Source/Site Ground Motion Comparisons here]("+bbpDir.getName()+"/"+sourceSiteLink+"/)");
				}
				if (!sourceDetailLinks.isEmpty()) {
					for (int i=0; i<sourceDetailLinks.size(); i++) {
						lines.add("");
						lines.add("["+sourceDetailNames.get(i)+"]("+bbpDir.getName()+"/"+sourceDetailLinks.get(i)+"/)");
					}
				}
			}
			if (!siteHazardLinksTable.isEmpty()) {
				lines.add("");
				lines.add("### Site Hazard Comparisons");
				lines.add(topLink);
				lines.add("");
				for (String gmpe : siteHazardLinksTable.keySet()) {
					List<String> siteNames = siteHazardNamesTable.get(gmpe);
					List<String> siteLinks = siteHazardLinksTable.get(gmpe);
					
					for (int i=0; i<siteNames.size(); i++)
						lines.add("* ["+siteNames.get(i)+"]("+siteLinks.get(i)+"/)");
				}
			}
			if (!partBLinks.isEmpty()) {
				lines.add("");
				lines.add("### BBP Part B Analysis, "+vm);
				lines.add(topLink);
				lines.add("");
				for (int i=0; i<partBLinks.size(); i++)
					lines.add("[BBP Part B Analysis Here, "+partBNames.get(i)+"]("+bbpDir.getName()+"/"+partBLinks.get(i)+")");
			}
			if (!rotatedRupLinks.isEmpty()) {
				lines.add("");
				lines.add("### Rotated Rupture Variability Comparisons, "+vm);
				lines.add(topLink);
				lines.add("");
				for (int i=0; i<rotatedRupLinks.size(); i++)
					lines.add("* ["+rotatedRupNames.get(i)+"]("+bbpDir.getName()+"/"+rotatedRupLinks.get(i)+"/)");
			}
			if (!azimuthalRupLinks.isEmpty()) {
				lines.add("");
				lines.add("### Scenario Spatial Distribution Plots, "+vm);
				lines.add(topLink);
				lines.add("");
				for (int i=0; i<azimuthalRupLinks.size(); i++)
					lines.add("* ["+azimuthalRupNames.get(i)+"]("+bbpDir.getName()+"/"+azimuthalRupLinks.get(i)+"/)");
			}
		}
		
		if (partBSummaryLink != null) {
			lines.add("");
			lines.add("## BBP PartB Summary");
			lines.add(topLink);
			lines.add("");
			lines.add("[BBP PartB Summary Here]("+partBSummaryLink+"/)");
		}
		
		if (crustalSubductionGMsLink != null) {
			lines.add("");
			lines.add("## BBP Crustal + Subduction Combined GMs");
			lines.add(topLink);
			lines.add("");
			lines.add("[BBP Crustal + Subduction Combined GMs Here]("+crustalSubductionGMsLink+"/)");
		}
		
		if (vmCompareRotRupLink != null) {
			lines.add("");
			lines.add("## BBP Velocity Model Comparisons");
			lines.add(topLink);
			lines.add("");
			lines.add("[BBP Velocity Model Comparisons Here]("+vmCompareRotRupLink+"/)");
		}
		
		if (plots) {
			File resourcesDir = new File(dir, "resources");
			Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
			lines.add("");
			int skipYears;
			double duration = getDurationYears();
			if (duration > 20000)
				skipYears = 5000;
			else if (duration > 10000)
				skipYears = 3000;
			else if (duration > 1000)
				skipYears = 1000;
			else
				skipYears = 0;
			lines.addAll(writeStandardDiagnosticPlots(resourcesDir, skipYears, plotMinMag, replot, topLink, standardPlots));
		}
		
		File inputFile = getParamFile();
		if (params != null) {
			lines.add("");
			lines.add("## Input File");
			lines.add(topLink);
			lines.add("");
			lines.add("```");
			for (String line : Files.readLines(inputFile, Charset.defaultCharset()))
				lines.add(line);
			lines.add("```");
		}
		
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		
		MarkdownUtils.writeReadmeAndHTML(lines, dir);
		
		// write metadata
		XMLUtils.writeObjectToXMLAsRoot(this, new File(dir, "catalog.xml"));
	}
	
	private String checkMoveFilter(File subDir, FilterMethod filter) throws IOException {
		String name = subDir.getName();
		if (!name.contains(filter.getPrefix())) {
			String target = name+"_filter_"+filter.getPrefix();
			System.out.println("Renaming (missing filter): "+name+" to "+target);
			Files.move(subDir, new File(subDir.getParentFile(), target));
			name = target;
		}
		return name;
	}

	public File getGeomFile() throws FileNotFoundException {
		// first try from the input file
		try {
			File file = new File(dir, getParams().get("faultFname"));
			if (file.exists())
				return file;
			System.out.println("Geometry file defined in input file (faultFname = "+file.getName()+
					") not found, searching for alternatives");
		} catch (Exception e) {
			System.out.println("Error locating geometry file name from parameter file: "+e.getMessage());
		}
		for (File file : dir.listFiles()) {
			String name = file.getName().toLowerCase();
			if (name.endsWith(".flt"))
				return file;
			if (name.startsWith("zfault") && name.endsWith(".in") && !name.contains("deepen_"))
				return file;
		}
		throw new FileNotFoundException("No geometry file found in "+dir.getAbsolutePath());
	}
	
	public synchronized List<SimulatorElement> getElements() throws IOException {
		if (elements == null) {
			File geomFile = getGeomFile();
			elements = RSQSimFileReader.readGeometryFile(geomFile, utmZone, utmBand);
		}
		return elements;
	}
	
	public Loader loader() throws IOException {
		return new Loader(getElements(), getCatalogDir());
	}
	
	private static File getTransFile(File dir, boolean transV) throws FileNotFoundException {
		for (File file : dir.listFiles()) {
			String name = file.getName().toLowerCase();
			if ((transV && name.startsWith("transv.") || !transV && name.startsWith("trans.")) && name.endsWith(".out"))
				return file;
		}
		throw new FileNotFoundException("No transitions file found in "+dir.getAbsolutePath()+". TransV? "+transV);
	}
	
	public synchronized RSQSimStateTransitionFileReader getTransitions() throws IOException {
		if (transReader == null) {
//			File transFile = getTransFile(getCatalogDir(), isVariableSlipSpeed());
			File transFile = null;
			TransVersion transVersion = null;
			boolean transV = isVariableSlipSpeed();
			Map<String, String> params = getParams();
			
			String newStyleStr = params.get("writeSlipSpeedInMainTransFile");
			if (newStyleStr != null) {
				int newStyleInt = Integer.parseInt(newStyleStr);
				if (newStyleInt == 0)
					transVersion = TransVersion.ORIGINAL;
				else
					transVersion = TransVersion.CONSOLIDATED_RELATIVE;
			}
			for (File file : getCatalogDir().listFiles()) {
				String name = file.getName().toLowerCase();
				if (!name.startsWith("trans") || !name.endsWith(".out"))
					continue;
				if (name.startsWith("transv.")) {
					transFile = file;
					Preconditions.checkState(transV, "have a transV file params indicate not transv: %s", file.getName());
					transVersion = TransVersion.TRANSV;
					// if we have a transv, use that. might have an old style plain trans as well
					break;
				} else {
					transFile = file;
				}
			}
//			System.out.println("Trans file: "+transFile.getAbsolutePath());
//			System.out.println("Trans version: "+transVersion);
//			System.out.println("TransV: "+transV);
			transReader = new RSQSimStateTransitionFileReader(transFile, getElements(), transVersion);
			transVersion = transReader.getVersion();
			if (transVersion == TransVersion.ORIGINAL)
				transReader.setPatchFixedVelocities(getSlipVelocities());
		}
		return transReader;
	}
	
	public synchronized RSQSimEventSlipTimeFunc getSlipTimeFunc(RSQSimEvent event) throws IOException {
		return new RSQSimEventSlipTimeFunc(getTransitions().getTransitions(event));
	}

	static GregorianCalendar cal(int year, int month, int day) {
		return new GregorianCalendar(year, month-1, day);
	}
	
	public double getMinSubSectFractForInclusion() {
		return minFractForInclusion;
	}
	
	public void setFractForInclusion(double minFractForInclusion) {
		if (subSectMapper != null)
			subSectMapper.setMinFractForInclusion(minFractForInclusion);
		this.minFractForInclusion = minFractForInclusion;
	}
	
	public synchronized List<? extends FaultSection> getSubSects() throws IOException {
		if (subSects == null)
			subSects = getDeformationModel().build(getFaultModel());
		return subSects;
	}
	
	private File getSolCacheDir() {
		File scratchDir = UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR;
		if (scratchDir.exists()) {
			// eclipse project
			File dir = new File(scratchDir, "ucerf3_fm_dm_sols");
			if (!dir.exists())
				Preconditions.checkState(dir.mkdir());
			return dir;
		} else {
			// use home dir
			String path = System.getProperty("user.home");
			File homeDir = new File(path);
			Preconditions.checkState(homeDir.exists(), "user.home dir doesn't exist: "+path);
			File openSHADir = new File(homeDir, ".opensha");
			if (!openSHADir.exists())
				Preconditions.checkState(openSHADir.mkdir(),
						"Couldn't create OpenSHA store location: "+openSHADir.getAbsolutePath());
			File uc3Dir = new File(openSHADir, "ucerf3_fm_dm_sols");
			if (!uc3Dir.exists())
				Preconditions.checkState(uc3Dir.mkdir(),
						"Couldn't create UCERF3 ERF store location: "+uc3Dir.getAbsolutePath());
			return uc3Dir;
		}
	}
	
	public FaultSystemSolution getCompareSol() throws IOException {
		synchronized (compSolsTable) {
			FaultSystemSolution sol = compSolsTable.get(fm, dm);
			
			if (sol == null) {
				File solDir = getSolCacheDir();
				File solFile = new File(solDir, fm.getFilePrefix()+"_"+dm.getFilePrefix()+"_MEAN_BRANCH_AVG_SOL.zip");
				if (!solFile.exists()) {
					// download it
					String addr = "http://opensha.usc.edu/ftp/kmilner/ucerf3/2013_05_10-ucerf3p3-production-10runs_fm_dm_sub_plots/"
							+ fm.getFilePrefix()+"_"+dm.getFilePrefix()+"/"+solFile.getName();
					FileUtils.downloadURL(addr, solFile);
				}
				
				sol = FaultSystemSolution.load(solFile);
				compSolsTable.put(fm, dm, sol);
			}
			
			return sol;
		}
	}
	
	private synchronized Map<Integer, Double> getSubSectAreas() throws IOException {
		if (subSectAreas == null)
			subSectAreas = RSQSimUtils.calcSubSectAreas(getElements(), getSubSects());
		return subSectAreas;
	}
	
	public synchronized Map<IDPairing, Double> getSubSectDistsCache() {
		if (subSectDistsCache == null)
			subSectDistsCache = new HashMap<>();
		return subSectDistsCache;
	}
	
	public synchronized RSQSimSubSectionMapper getSubSectMapper() throws IOException {
		if (subSectMapper == null)
			subSectMapper = new RSQSimSubSectionMapper(getSubSects(), getElements(), minFractForInclusion,
					getSubSectAreas(), getSubSectDistsCache());
		return subSectMapper;
	}
	
	public RSQSimEqkRupture getEqkRupture(RSQSimEvent event) {
//		if (getFaultModel() != null) {
//			try {
//				getMappedSubSectRupture(event);
//			} catch (Exception e) {}
//		}
		return getCumDistCalcRupture(event);
	}
	
	private LocationElementDistanceCacheFactory locElemDistCacheFactory = null;
	
	public RSQSimEqkRupture getCumDistCalcRupture(RSQSimEvent event) {
		if (locElemDistCacheFactory == null) {
			synchronized (this) {
				if (locElemDistCacheFactory == null)
					locElemDistCacheFactory = new LocationElementDistanceCacheFactory();
			}
		}
		return RSQSimUtils.buildCumDistRupture(event, locElemDistCacheFactory);
	}
	
	public RSQSimSubSectEqkRupture getMappedSubSectRupture(RSQSimEvent event) {
		RSQSimSubSectionMapper mapper;
		try {
			mapper = getSubSectMapper();
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		return RSQSimUtils.buildSubSectBasedRupture(mapper, event);
	}
	
	public List<FaultSection> getSubSectsForRupture(RSQSimEvent event) {
		RSQSimSubSectionMapper mapper;
		try {
			mapper = getSubSectMapper();
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		List<List<SubSectionMapping>> bundled =  mapper.getFilteredSubSectionMappings(event);
		if (minFractForInclusion >= 0 && bundled.isEmpty())
			bundled = mapper.getAllSubSectionMappings(event);
		List<FaultSection> allSects = new ArrayList<>();
		for (List<SubSectionMapping> bundle : bundled)
			for (SubSectionMapping mapping : bundle)
				allSects.add(mapping.getSubSect());
		return allSects;
	}
	
	public class Loader {
		private List<SimulatorElement> elements;
		private File catalogDir;
		private boolean skipSlipsAndTimes = false;
		
		private List<RuptureIdentifier> loadIdens;
		
		private Loader(List<SimulatorElement> elements, File catalogDir) {
			super();
			this.elements = elements;
			this.catalogDir = catalogDir;
			
			loadIdens = new ArrayList<>();
		}
		
		public Loader magRange(double minMag, double maxMag) {
			loadIdens.add(new MagRangeRuptureIdentifier(minMag, maxMag));
			return this;
		}
		
		public Loader minMag(double minMag) {
			return magRange(minMag, Double.POSITIVE_INFINITY);
		}
		
		public Loader maxMag(double maxMag) {
			return magRange(Double.NEGATIVE_INFINITY, maxMag);
		}
		
		public Loader skipSlipsAndTimes() {
			this.skipSlipsAndTimes = true;
			return this;
		}
		
		public Loader skipYears(double years) {
			if (years > 0)
				loadIdens.add(new SkipYearsLoadIden(years));
			return this;
		}
		
		public Loader withinTimeRange(double tStartSecs, double tEndSecs) {
			loadIdens.add(new EventTimeIdentifier(tStartSecs, tEndSecs, false));
			return this;
		}
		
		public Loader maxDuration(double years) {
			loadIdens.add(new CatalogLengthLoadIden(years));
			return this;
		}
		
		public Loader matches(RuptureIdentifier iden) {
			loadIdens.add(iden);
			return this;
		}
		
		public Loader withinCutoffDist(double maxDist, Collection<Location> locs) {
			return withinCutoffDist(maxDist, locs.toArray(new Location[0]));
		}
		
		public Loader withinCutoffDist(double maxDist, Location... locs) {
			List<RegionIden> regIdens = new ArrayList<>();
			for (Location loc : locs)
				regIdens.add(new RegionIden(new Region(loc, maxDist)));
			return matches(new LogicalOrRupIden(regIdens));
		}
		
		public Loader hasTransitions() throws IOException {
			loadIdens.add(new RSQSimTransValidIden(getTransitions()));
			return this;
		}
		
		public Loader forParentSections(boolean calcU3Offset, int... parentIDs) throws IOException {
			List<Integer> sectionIDs = new ArrayList<>();
			for (FaultSection sect : getSubSects())
				if (Ints.contains(parentIDs, sect.getParentSectionId()))
					sectionIDs.add(sect.getSectionId());
			return forSections(calcU3Offset, Ints.toArray(sectionIDs));
		}
		
		public Loader forSections(boolean calcU3Offset, int... sectionIDs) throws IOException {
			if (calcU3Offset) {
				int offset = RSQSimUtils.getSubSectIndexOffset(getElements(), getSubSects());
				if (offset != 0) {
					for (int i=0; i<sectionIDs.length; i++)
						sectionIDs[i] -= offset;
				}
			}
			loadIdens.add(new SectionIDIden(XML_METADATA_NAME, getElements(), sectionIDs));
			return this;
		}
		
		public RSQSimEvent byID(int eventID) throws IOException {
			List<RSQSimEvent> events = this.byIDs(eventID);
			Preconditions.checkState(events.size() == 1, "Event "+eventID+" not found");
			return events.get(0);
		}
		
		public List<RSQSimEvent> byIDs(int... eventIDs) throws IOException {
			loadIdens.add(new EventIDsRupIden(eventIDs));
			return this.load();
		}
		
		public Loader maxEventID(int maxID) {
			loadIdens.add(new EventIDsRangeIden(com.google.common.collect.Range.closed(Integer.MIN_VALUE, maxID)));
			return this;
		}
		
		public Loader minEventID(int minID) {
			loadIdens.add(new EventIDsRangeIden(com.google.common.collect.Range.closed(minID, Integer.MAX_VALUE)));
			return this;
		}
		
		public Loader forEventIDRange(com.google.common.collect.Range<Integer> range) {
			loadIdens.add(new EventIDsRangeIden(range));
			return this;
		}
		
		public List<RSQSimEvent> load() throws IOException {
			LogicalAndRupIden loadIden = new LogicalAndRupIden(loadIdens);
			List<RuptureIdentifier> rupIdens = new ArrayList<>();
			rupIdens.add(loadIden);
			return RSQSimFileReader.readEventsFile(catalogDir, elements, rupIdens, skipSlipsAndTimes);
		}
		
		public Iterable<RSQSimEvent> iterable() throws IOException {
			LogicalAndRupIden loadIden = new LogicalAndRupIden(loadIdens);
			List<RuptureIdentifier> rupIdens = new ArrayList<>();
			rupIdens.add(loadIden);
			return RSQSimFileReader.getEventsIterable(catalogDir, elements, rupIdens, skipSlipsAndTimes);
		}
		
		public Map<Integer, RSQSimEvent> loadMap() throws IOException {
			List<RSQSimEvent> events = load();
			Map<Integer, RSQSimEvent> map = new HashMap<>(events.size());
			for (RSQSimEvent event : events)
				map.put(event.getID(), event);
			return map;
		}
	}
	
	public FaultSystemSolution buildSolution(Loader loader, double minMag) throws IOException {
		return buildSolution(loader.load(), minMag);
	}
	
	public FaultSystemSolution buildSolution(List<RSQSimEvent> events, double minMag) throws IOException {
		return RSQSimUtils.buildFaultSystemSolution(getSubSects(), getElements(), events, minMag, minFractForInclusion);
	}
	
	// TODO duplicate with getCompareSol()?
	public FaultSystemSolution getComparisonSolution() throws IOException {
		File solFile = new File(fmDmSolDir, getFaultModel().getFilePrefix()
				+"_"+getDeformationModel().getFilePrefix()+"_MEAN_BRANCH_AVG_SOL.zip");
		if (solFile.exists()) {
			System.out.println("Loading comparison FSS from "+solFile.getAbsolutePath());
			compSol = FaultSystemSolution.load(solFile);
		} else {
			System.out.println("Comparison sol file doesn't exist: "+solFile.getAbsolutePath());
		}
		return compSol;
	}
	
	public enum StandardPlots {
		MFD,
		MAG_AREA,
		SLIP_AREA,
		SLIP_LEN,
		SLIP_ALONG_RUPTURE,
		SLIP_RATE,
		RUP_VELOCITY,
		GLOBAL_RECURRENCE,
		NORM_RECURRENCE,
		U3_NORM_RECURRENCE,
		STATIONARITY,
		SECTION_RECURRENCE,
		SECT_PARTIC_NUCL,
		PALEO_RECURRENCE,
		PALEO_OPEN_INTERVAL,
		MOMENT_RATE_VARIABILITY,
		ELASTIC_REBOUND_TRIGGERING
	}
	
	private List<AbstractPlot> buildStandardPlotLines(File outputDir, int skipYears, double minMag, boolean replot, String topLink, List<String> lines,
			HashSet<StandardPlots> plotsSet) throws IOException {
		lines.add("## Plots");
		List<AbstractPlot> plots = new ArrayList<>();
		
		boolean hasFM = getFaultModel() != null && getDeformationModel() != null;
		
		TableBuilder table;
		
		if (plotsSet.contains(StandardPlots.MFD)) {
			if (replot || !new File(outputDir, "mfd.png").exists()) {
				MFDPlot mfdPlot = new MFDPlot(minMag);
				File mfdCSV = null;
				if (hasFM)
					mfdCSV = new File(fmDmSolDir, getFaultModel().getFilePrefix()
							+"_"+getDeformationModel().getFilePrefix()+"_supra_plus_sub_seis_cumulative.csv");
//				System.out.println(mfdCSV.getAbsolutePath()+" ? "+mfdCSV.exists());
				if (mfdCSV != null && mfdCSV.exists()) {
					System.out.println("Loading UCERF3 comparison MFD: "+mfdCSV.getAbsolutePath());
					CSVFile<String> csv = CSVFile.readFile(mfdCSV, true);
					DiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
					DiscretizedFunc minFunc = new ArbitrarilyDiscretizedFunc();
					DiscretizedFunc maxFunc = new ArbitrarilyDiscretizedFunc();
					for (int row=1; row<csv.getNumRows(); row++) {
						double mag = Double.parseDouble(csv.get(row, 0));
						if (mag < minMag-0.01)
							continue;
						double mean = Double.parseDouble(csv.get(row, 1));
						if (mean == 0d)
							break;
						meanFunc.set(mag, mean);
						minFunc.set(mag, Double.parseDouble(csv.get(row, 2)));
						maxFunc.set(mag, Double.parseDouble(csv.get(row, 3)));
					}
					UncertainArbDiscFunc compRange = new UncertainArbDiscFunc(meanFunc, minFunc, maxFunc);
					compRange.setName("U3 On Fault");
					mfdPlot.setComparableRange(compRange);
				}
				mfdPlot.initialize(getName(), outputDir, "mfd");
				plots.add(mfdPlot);
			}
			lines.add("### Magnitude-Frequency Plot");
			lines.add(topLink);
			lines.add("");
			lines.add("![MFD]("+outputDir.getName()+"/mfd.png)");
		}
		
		if (plotsSet.contains(StandardPlots.MAG_AREA)) {
			if (replot || !new File(outputDir, "mag_area_hist2D.png").exists()) {
				MagAreaScalingPlot magAreaPlot = new MagAreaScalingPlot(false);
				magAreaPlot.initialize(getName(), outputDir, "mag_area");
				plots.add(magAreaPlot);
			}
			lines.add("### Magnitude-Area Plots");
			lines.add(topLink);
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.addLine("Scatter", "2-D Hist");
			table.initNewLine();
			table.addColumn("![MA Scatter]("+outputDir.getName()+"/mag_area.png)");
			table.addColumn("![MA Hist]("+outputDir.getName()+"/mag_area_hist2D.png)");
			table.finalizeLine();
			lines.addAll(table.build());
			
			// look for fault style specific plots
			table = null;
			for (FaultStyle style : FaultStyle.values()) {
				File scatterPlot = new File(outputDir, "mag_area_"+style.name()+".png");
				if (!scatterPlot.exists())
					continue;
				if (table == null) {
					table = MarkdownUtils.tableBuilder();
					table.addLine("Fault Style", "Scatter", "2-D Hist");
				}
				table.initNewLine();
				table.addColumn("**"+style+"**");
				table.addColumn("![MA Scatter]("+outputDir.getName()+"/mag_area_"+style.name()+".png)");
				table.addColumn("![MA Hist]("+outputDir.getName()+"/mag_area_"+style.name()+"_hist2D.png)");
				table.finalizeLine();
			}
			if (table != null) {
				lines.add("#### Mechanism-Dependent Magnitude-Area Plots");
				lines.add(topLink);
				lines.add("");
				lines.add("Here we disaggregate the magnitude-area scaling plots by focal mechanism. Multi-fault ruptures which incorporate "
						+ "multiple faulting styles are included in plot for the dominent fault style so long as no more than 10% of the "
						+ "participating elements are of a different style, otherwise they are listed as 'Unknown'.");
				lines.add("");
				lines.addAll(table.build());
			}
		}
		
		if (plotsSet.contains(StandardPlots.SLIP_AREA)) {
			if (replot || !new File(outputDir, "slip_area_hist2D.png").exists()) {
				MagAreaScalingPlot magAreaPlot = new MagAreaScalingPlot(true);
				magAreaPlot.initialize(getName(), outputDir, "slip_area");
				plots.add(magAreaPlot);
			}
			lines.add("### Slip-Area Plots");
			lines.add(topLink);
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.addLine("Scatter", "2-D Hist");
			table.initNewLine();
			table.addColumn("![Slip Area Scatter]("+outputDir.getName()+"/slip_area.png)");
			table.addColumn("![Slip Area Hist]("+outputDir.getName()+"/slip_area_hist2D.png)");
			table.finalizeLine();
			lines.addAll(table.build());
		}

		SlipAlongSectAlgorithm defaultSlipAlg = SlipAlongSectAlgorithm.MID_SEIS_SLIPPED_LEN;
		String midSeisDescription = "no deeper than "+optionalDigitDF.format(RSQSimSubSectionMapper.MID_SEIS_MAX_DEPTH_DEFAULT)+" km, "
				+ "no shallower than "+optionalDigitDF.format(RSQSimSubSectionMapper.MID_SEIS_MIN_DEPTH_DEFAULT)+" km, "
				+ "and no less than "+optionalDigitDF.format(RSQSimSubSectionMapper.MID_SEIS_BUFFER_DEFAULT)+" km down- or "
				+ "up-dip from the top or bottom of the fault";
		
		if (plotsSet.contains(StandardPlots.SLIP_LEN) && hasFM) {
			if (replot || !new File(outputDir, "slip_len_"+SlipAlongSectAlgorithm.MID_SEIS_SURF_SLIP_LEN.name()+".png").exists()) {
				SlipLengthScalingPlot slipLengthPlot = new SlipLengthScalingPlot(getSubSectMapper(), 6.5);
				slipLengthPlot.setDisaggregateFaultStyles(defaultSlipAlg);
				slipLengthPlot.initialize(getName(), outputDir, "slip_len");
				plots.add(slipLengthPlot);
			}
			lines.add("### Slip-Length Plots");
			lines.add(topLink);
			lines.add("");
			lines.add("These plots compute average slip-length scaling at mid-seismogenic depth. We define mid-seismogenic depth " 
					+ "to be "+midSeisDescription+". Average slip is computed across all elements in this "
					+ "mid-seismogenic region, including any which did not slip, along the full length of the rupture.");
			lines.add("");
			lines.add("We define the rupture length, which also determines the region at mid-seismogenic depth across which we "
					+ "compute average slip, multiple ways in order to test sensitivity:");
			lines.add("");
			for (SlipAlongSectAlgorithm slipAlg : SlipAlongSectAlgorithm.values())
				lines.add("* **"+slipAlg+":** "+slipAlg.getDescription());
			lines.add("");
			File examplePlot = new File(outputDir, "slip_len_example_rupture.png");
			if (examplePlot.exists()) {
				lines.add("These length algorithms are illustrated in the following example plot, which also has the "
						+ "mid-seismogenic depth range outlined in a cyan dashed line:");
				lines.add("");
				lines.add("![Example plot]("+outputDir.getName()+"/"+examplePlot.getName()+")");
				lines.add("");
			}
			lines.add("The average value is plotted in a thick gray line, and UCERF3 Scaling Relationships in colored lines "
					+ "(assuming a down dip width of "+optionalDigitDF.format(SlipLengthScalingPlot.LEN_COMP_DDW)+" km).");
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.addLine("Lengh Algorithm", "Scatter", "2-D Hist");
			for (SlipAlongSectAlgorithm slipAlg : SlipAlongSectAlgorithm.values()) {
				table.initNewLine();
				table.addColumn("**"+slipAlg+"**");
				String prefix = "slip_len_"+slipAlg.name();
				table.addColumn("![Slip Length Scatter]("+outputDir.getName()+"/"+prefix+".png)");
				table.addColumn("![Slip Length Hist]("+outputDir.getName()+"/"+prefix+"_hist2D.png)");
				table.finalizeLine();
			}
			lines.addAll(table.build());
			
			// look for fault style specific plots
			table = null;
			for (FaultStyle style : FaultStyle.values()) {
				File scatterPlot = new File(outputDir, "mag_area_"+style.name()+".png");
				if (!scatterPlot.exists())
					continue;
				if (table == null) {
					table = MarkdownUtils.tableBuilder();
					table.addLine("Fault Style", "Scatter", "2-D Hist");
				}
				table.initNewLine();
				table.addColumn("**"+style+"**");
				String prefix = "slip_len_"+defaultSlipAlg.name()+"_"+style.name();
				table.addColumn("![Slip Length Scatter]("+outputDir.getName()+"/"+prefix+".png)");
				table.addColumn("![Slip Length Hist]("+outputDir.getName()+"/"+prefix+"_hist2D.png)");
				table.finalizeLine();
			}
			if (table != null) {
				lines.add("#### Mechanism-Dependent Slip-Length Plots");
				lines.add(topLink);
				lines.add("");
				lines.add("Here we disaggregate the slip-length scaling plots by focal mechanism. Multi-fault ruptures which incorporate "
						+ "multiple faulting styles are included in plot for the dominent fault style so long as no more than 10% of the "
						+ "participating elements are of a different style, otherwise they are listed as 'Unknown'.");
				lines.add("");
				lines.addAll(table.build());
			}
		}
		
		List<Range> lengthBins = new ArrayList<>();
		lengthBins.add(null);
		lengthBins.add(new Range(0d, 25));
		lengthBins.add(new Range(25, 50));
		lengthBins.add(new Range(50, 100));
		lengthBins.add(new Range(100, Double.POSITIVE_INFINITY));
		if (plotsSet.contains(StandardPlots.SLIP_ALONG_RUPTURE) && hasFM) {
			if (replot || !new File(outputDir, "slip_along_rupture_multi_norm.png").exists()) {
				SlipAlongRupturePlot slipAlongPlot = new SlipAlongRupturePlot(getSubSectMapper(), 6.5, defaultSlipAlg, fm.getNamedFaults(), lengthBins);
				slipAlongPlot.initialize(getName(), outputDir, "slip_along_rupture");
				plots.add(slipAlongPlot);
			}
			lines.add("### Slip Along Rupture (Dsr) Plots");
			lines.add(topLink);
			lines.add("");
			lines.add("These plots show the slip along rupture distiribution, noted D<sub>SR</sub> in UCERF3. First we compute average "
					+ "slip along each mapped subsection at mid-seismogenic depth (using the *"+defaultSlipAlg+"* algorithm), then plot that slip along "
					+ "strike, normalized by the average slip across all subsections in that rupture. We do this for single-fault events, which "
					+ "can span multiple segments (e.g. SAF Mojave and San Bernardino), and also separately for each junction in multi-fault events. "
					+ "This is done using the UCERF3 'named faults' list to determine if multiple fault sections belong to the same master fault. "
					+ "We only consider ruptures where at least 2 subsections participated (2 on each side of the jump for multi-fault ruptures).");
			lines.add("");
			lines.add("Ruptures are binned by their length in each row below. For multi-fault ruptures, the junction point is at x=0 with the shorter "
					+ "side of the rupture on the left (below zero), and longer half on the right");
			lines.add("");
			lines.add("Average values are plotted with a solid black line, and "+(float)+SlipAlongRupturePlot.SQRT_SINE_SCALAR+"*sqrt(sin(|x*&pi;|)) "
					+ "in a dashed gray line (normalized length plots only).");
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.addLine("Rupture Length", "Single-fault, absolute distance from either rupture endpoint",
					"Single-fault, normalized distance along strike", "Multi-fault, normalized distance on either side of jump");
			for (Range lengthBin : lengthBins) {
				table.initNewLine();
				String lenStr;
				String prefixAdd = "";
				if (lengthBin == null) {
					lenStr = "All Lengths";
				} else {
					if (Double.isInfinite(lengthBin.getUpperBound())) {
						lenStr = "Len≥"+optionalDigitDF.format(lengthBin.getLowerBound());
						prefixAdd = "_len_"+optionalDigitDF.format(lengthBin.getLowerBound())+"+";
					} else {
						lenStr = "Len=["+optionalDigitDF.format(lengthBin.getLowerBound())+" "+optionalDigitDF.format(lengthBin.getUpperBound())+"]";
						prefixAdd = "_len_"+optionalDigitDF.format(lengthBin.getLowerBound())+"_"+optionalDigitDF.format(lengthBin.getUpperBound());
					}
				}
				table.addColumn("**"+lenStr+"**");
				table.addColumn(imageIfExists(outputDir, "slip_along_rupture_single_abs"+prefixAdd+".png", "Slip Along Rupture", "N/A"));
				table.addColumn(imageIfExists(outputDir, "slip_along_rupture_single_norm"+prefixAdd+".png", "Slip Along Rupture", "N/A"));
				table.addColumn(imageIfExists(outputDir, "slip_along_rupture_multi_norm"+prefixAdd+".png", "Slip Along Rupture", "N/A"));
				table.finalizeLine();
			}
			lines.addAll(table.build());
			lines.add("");
			lines.add("#### Two- and Three-Fault Slip Along Rupture");
			lines.add("");
			lines.add("These plots show D<sub>SR</sub> for two- and three-fault ruptures. Lengths are normalized, with the first fault in x=[0 1], "
					+ "second in x=[1 2], etc. Rupture are organized such that the leftmost side is always shorter than the rightmost side.");
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.addLine("Rupture Length", "Two-fault Ruptures", "Three-Fault Ruptures");
			for (Range lengthBin : lengthBins) {
				table.initNewLine();
				String lenStr;
				String prefixAdd = "";
				if (lengthBin == null) {
					lenStr = "All Lengths";
				} else {
					if (Double.isInfinite(lengthBin.getUpperBound())) {
						lenStr = "Len≥"+optionalDigitDF.format(lengthBin.getLowerBound());
						prefixAdd = "_len_"+optionalDigitDF.format(lengthBin.getLowerBound())+"+";
					} else {
						lenStr = "Len=["+optionalDigitDF.format(lengthBin.getLowerBound())+" "+optionalDigitDF.format(lengthBin.getUpperBound())+"]";
						prefixAdd = "_len_"+optionalDigitDF.format(lengthBin.getLowerBound())+"_"+optionalDigitDF.format(lengthBin.getUpperBound());
					}
				}
				table.addColumn("**"+lenStr+"**");
				table.addColumn(imageIfExists(outputDir, "slip_along_rupture_two_norm"+prefixAdd+".png", "Slip Along Rupture", "N/A"));
				table.addColumn(imageIfExists(outputDir, "slip_along_rupture_three_norm"+prefixAdd+".png", "Slip Along Rupture", "N/A"));
				table.finalizeLine();
			}
			lines.addAll(table.build());
		}
		
		if (plotsSet.contains(StandardPlots.SLIP_RATE) && hasFM) {
			if (replot || !new File(outputDir, "slip_rate_table.csv").exists()) {
				SlipRateComparePlot plot = new SlipRateComparePlot(getSubSectMapper(), getFaultModel(), getDeformationModel(),
						getComparisonSolution());
				plot.initialize(getName(), outputDir, "slip_rate");
				plots.add(plot);
			}
			
			lines.add("### Slip Rate Plots");
			lines.add(topLink);
			lines.add("");
			String slipRateDesc = "Slip rates are calculated at mid-seismogenic depth: "+midSeisDescription+". UCERF3 comparisons are included "
					+ "with the original target slip rate for the fault and deformation model used as input to the simulator when constructing "
					+ "the geometry, but this target is often smoothed and/or modified before use in the simlators.";
			boolean hasSlipU3Sol = new File(outputDir, "slip_rate_u3_sol_map.png").exists();
			if (hasSlipU3Sol)
				slipRateDesc += " Post-UCERF3 inversion slip rates (which will not perfectly match the target) are also included and labeled as "
					+ "'UCERF3 Solution'.";
			lines.add(slipRateDesc);
			lines.add("");
			String csvLink = outputDir.getName()+"/slip_rate_table.csv";
			lines.add("Download slip rate data: ["+csvLink+"]("+csvLink+")");
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			table.addColumn("<p align=\"center\">**Simulation Slip Rate**</p>");
			table.addColumn("<p align=\"center\">**Simulation vs Target Ratio**</p>");
			table.addColumn("<p align=\"center\">**UCERF3 Target Slip Rate**</p>");
			table.addColumn("<p align=\"center\">**Simulation vs UCERF3 Target Ratio**</p>");
			if (hasSlipU3Sol) {
				table.addColumn("<p align=\"center\">**UCERF3 Solution Slip Rate**</p>");
				table.addColumn("<p align=\"center\">**UCERF3 Solution vs Target Ratio**</p>");
			}
			table.finalizeLine();
			table.initNewLine();
			table.addColumn("![Slip Rate Plot]("+outputDir.getName()+"/slip_rate_sim_map.png)");
			table.addColumn("![Slip Rate Plot]("+outputDir.getName()+"/slip_rate_sim_ratio_map.png)");
			table.addColumn("![Slip Rate Plot]("+outputDir.getName()+"/slip_rate_u3_target_map.png)");
			table.addColumn("![Slip Rate Plot]("+outputDir.getName()+"/slip_rate_sim_u3_ratio_map.png)");
			if (hasSlipU3Sol) {
				table.addColumn("![Slip Rate Plot]("+outputDir.getName()+"/slip_rate_u3_sol_map.png)");
				table.addColumn("![Slip Rate Plot]("+outputDir.getName()+"/slip_rate_u3_ratio_map.png)");
			}
			table.finalizeLine();
			lines.addAll(table.wrap(4, 0).build());
			
			NamedFaults named = getFaultModel().getNamedFaults();
			if (named != null) {
				Map<String, String> slipFaultPlotMap = new HashMap<>();
				for (String fault : named.getFaultNames()) {
					File faultPlot = new File(outputDir, "slip_rate_fault_"+fault.replaceAll("\\W+", "_")+".png");
					if (faultPlot.exists())
						slipFaultPlotMap.put(fault, faultPlot.getName());
				}
				if (slipFaultPlotMap.size() > 0) {
					lines.add("#### Slip Rate Fault Plots");
					lines.add(topLink);
					lines.add("");
					table = MarkdownUtils.tableBuilder();
					List<String> sortedNames = ComparablePairing.getSortedData(slipFaultPlotMap);
					table.initNewLine();
					for (String faultName : sortedNames)
						table.addColumn("<p align=\"center\">**"+faultName+"**</p>");
					table.finalizeLine();
					table.initNewLine();
					for (String faultName : sortedNames)
						table.addColumn("![Slip Rate Plot]("+outputDir.getName()+"/"+slipFaultPlotMap.get(faultName)+")");
					table.finalizeLine();
					lines.addAll(table.wrap(3, 0).build());
				}
			}
		}
		
		if (plotsSet.contains(StandardPlots.RUP_VELOCITY)) {
			if (replot || !new File(outputDir, "rupture_velocity_scatter.png").exists()) {
				RuptureVelocityPlot rupVelPlot = new RuptureVelocityPlot(getElements(), minMag);
				rupVelPlot.initialize(getName(), outputDir, "rupture_velocity");
				plots.add(rupVelPlot);
			}
			lines.add("### Rupture Velocity Plots");
			lines.add(topLink);
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.initNewLine().addColumn("**Scatter**");
			table.addColumn("![Rupture Velocity Scatter]("+outputDir.getName()+"/rupture_velocity_scatter.png)");
			table.finalizeLine().initNewLine().addColumn("**Distance/Velocity**");
			table.addColumn("![Rupture Velocity vs Dist]("+outputDir.getName()+"/rupture_velocity_vs_dist.png)");
			table.finalizeLine();
			lines.addAll(table.build());
		}
		
		double[] riMinMags = {6d, 6.5, 7d, 7.5};
		while (minMag > riMinMags[0])
			riMinMags = Arrays.copyOfRange(riMinMags, 1, riMinMags.length);
		if (plotsSet.contains(StandardPlots.GLOBAL_RECURRENCE)) {
			if (replot || !new File(outputDir, "interevent_times_m7.5.png").exists()) {
				RecurrenceIntervalPlot riPlot = new RecurrenceIntervalPlot(riMinMags);
				riPlot.initialize(getName(), outputDir, "interevent_times");
				plots.add(riPlot);
			}
			lines.add("### Global Interevent-Time Distributions");
			lines.add(topLink);
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			for (double riMinMag : riMinMags)
				if (riMinMag == Math.round(riMinMag))
					table.addColumn("**M≥"+(int)riMinMag+"**");
				else
					table.addColumn("**M≥"+(float)riMinMag+"**");
			table.finalizeLine().initNewLine();
			for (double riMinMag : riMinMags)
				if (riMinMag == Math.round(riMinMag))
					table.addColumn("![Interevent Times]("+outputDir.getName()+"/interevent_times_m"+(int)riMinMag+".png)");
				else
					table.addColumn("![Interevent Times]("+outputDir.getName()+"/interevent_times_m"+(float)riMinMag+".png)");
			table.finalizeLine();
			lines.addAll(table.build());
		}
		
		if (plotsSet.contains(StandardPlots.NORM_RECURRENCE)) {
			if (replot || !new File(outputDir, "norm_ri_elem_m7.5.png").exists()) {
				List<NormalizedFaultRecurrenceIntervalPlot> myPlots = new ArrayList<>();
				myPlots.add(new NormalizedFaultRecurrenceIntervalPlot(getElements(), riMinMags));
				if (hasFM) {
					RSQSimSubSectionMapper mapper = getSubSectMapper();
					myPlots.add(new NormalizedFaultRecurrenceIntervalPlot(getElements(), SectType.SUBSECTION,
							mapper, riMinMags));
					myPlots.add(new NormalizedFaultRecurrenceIntervalPlot(getElements(), SectType.PARENT,
							mapper, riMinMags));
				}
				for (NormalizedFaultRecurrenceIntervalPlot plot : myPlots)
					plot.initialize(getName(), outputDir, "norm_ri_"+plot.getSectType().getPrefix());
				plots.addAll(myPlots);
			}
			lines.add("### Normalized Fault Interevent-Time Distributions");
			lines.add(topLink);
			lines.add("");
			lines.add("These plots show interevent-time distributions for a point on a fault (either an element, "
					+ " or aggregated at the subsection or parent section level).");
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			table.addColumn("");
			for (double riMinMag : riMinMags)
				if (riMinMag == Math.round(riMinMag))
					table.addColumn("**M≥"+(int)riMinMag+"**");
				else
					table.addColumn("**M≥"+(float)riMinMag+"**");
			table.finalizeLine();
			SectType[] types;
			if (hasFM)
				types = new SectType[] {SectType.ELEMENT, SectType.SUBSECTION, SectType.PARENT };
			else
				types = new SectType[] {SectType.ELEMENT};
			for (SectType type : types) {
				table.initNewLine();
				table.addColumn("**"+type.getSimType()+"s**");
				for (double riMinMag : riMinMags)
					if (riMinMag == Math.round(riMinMag))
						table.addColumn("![Norm RIs]("+outputDir.getName()+"/norm_ri_"+type.getPrefix()+"_m"+(int)riMinMag+".png)");
					else
						table.addColumn("![Norm RIs]("+outputDir.getName()+"/norm_ri_"+type.getPrefix()+"_m"+(float)riMinMag+".png)");
				table.finalizeLine();
			}
			lines.addAll(table.build());
		}

		BPTAveragingTypeOptions aveType = BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE;
		if (plotsSet.contains(StandardPlots.U3_NORM_RECURRENCE) && hasFM) {
			if (replot || !new File(outputDir, "u3_norm_ri_m7.5.png").exists()) {
				U3StyleNormalizedRuptureRecurrenceIntervalPlot plot = new U3StyleNormalizedRuptureRecurrenceIntervalPlot(
						getElements(), aveType, getSubSectMapper(), riMinMags);
				plot.initialize(getName(), outputDir, "u3_norm_ri");
				plots.add(plot);
			}
			lines.add("### Normalized Rupture Interevent-Time Distributions");
			lines.add(topLink);
			lines.add("");
			String line = "These plots show interevent-time distributions, averaged over a rupture, similar to "
					+ "the UCERF3 BPT calculation. For each rupture, we compute ";
			if (aveType.isAveNTS()) {
				line += "the average normalized open interval across all subsections which participate.";
			} else {
				line += "the average open interval across all subsections which participate, normalized "
						+ "by the average recurrence interval on those sections";
				if (aveType.isAveRI())
					line += ".";
				else
					line += ", which is computed as the inverse of the average rate on each subsection";
			}
			lines.add(line);
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			for (double riMinMag : riMinMags)
				if (riMinMag == Math.round(riMinMag))
					table.addColumn("**M≥"+(int)riMinMag+"**");
				else
					table.addColumn("**M≥"+(float)riMinMag+"**");
			table.finalizeLine();
			table.initNewLine();
			for (double riMinMag : riMinMags)
				if (riMinMag == Math.round(riMinMag))
					table.addColumn("![Norm RIs]("+outputDir.getName()+"/u3_norm_ri_m"+(int)riMinMag+".png)");
				else
					table.addColumn("![Norm RIs]("+outputDir.getName()+"/u3_norm_ri_m"+(float)riMinMag+".png)");
			table.finalizeLine();
			lines.addAll(table.build());
		}
		
		if (plotsSet.contains(StandardPlots.STATIONARITY)) {
			if (replot || !new File(outputDir, "stationarity.png").exists()) {
				StationarityPlot stationarityPlot = new StationarityPlot(minMag, 7d);
				stationarityPlot.initialize(getName(), outputDir, "stationarity");
				plots.add(stationarityPlot);
			}
			lines.add("### Stationarity Plot");
			lines.add(topLink);
			lines.add("");
			lines.add("![Stationarity]("+outputDir.getName()+"/stationarity.png)");
		}
		
		String testMagStr;
		if (riMinMags[0] == Math.floor(riMinMags[0]))
			testMagStr = (int)riMinMags[0]+"";
		else
			testMagStr = (float)riMinMags[0]+"";
		if (plotsSet.contains(StandardPlots.SECT_PARTIC_NUCL) && hasFM) {
			if (replot || !new File(outputDir, "sub_sect_m"+testMagStr+"_partic_map.png").exists()) {
				RSQSimSubSectionMapper mapper = getSubSectMapper();
				SectParticipationNucleationPlot sectPlot = new SectParticipationNucleationPlot(mapper, riMinMags);
				sectPlot.initialize(getName(), outputDir, "sub_sect");
				plots.add(sectPlot);
			}
			lines.add("### Subsection Participation & Nucleation Maps");
			lines.add(topLink);
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.addLine("Min Mag", "Participation Rate", "Nucleation Rate", "Nucleation/Participation Ratio");
			for (double riMinMag : riMinMags) {
				String prefix = "sub_sect_m";
				if (riMinMag == Math.floor(riMinMag))
					prefix += (int)riMinMag;
				else
					prefix += (float)riMinMag;
				table.initNewLine();
				table.addColumn("**M≥"+(float)riMinMag+"**");
				table.addColumn("![Partic]("+outputDir.getName()+"/"+prefix+"_partic_map.png)");
				table.addColumn("![Nucl]("+outputDir.getName()+"/"+prefix+"_nucl_map.png)");
				table.addColumn("![Ratio]("+outputDir.getName()+"/"+prefix+"_nucl_partic_ratio_map.png)");
				table.finalizeLine();
			}
			lines.addAll(table.build());
		}
		
		if (plotsSet.contains(StandardPlots.SECTION_RECURRENCE) && hasFM) {
			if (replot || !new File(outputDir, "interevent_elements_m"+testMagStr+"_scatter.png").exists()) {
				RSQSimSubSectionMapper mapper = getSubSectMapper();
				SectionRecurrenceComparePlot elemCompare = new SectionRecurrenceComparePlot(getElements(), getCompareSol(), "UCERF3",
						SectionRecurrenceComparePlot.SectType.ELEMENT, mapper, riMinMags);
				elemCompare.initialize(getName(), outputDir, "interevent_elements");
				plots.add(elemCompare);
				
				SectionRecurrenceComparePlot subSectCompare = new SectionRecurrenceComparePlot(getElements(), getCompareSol(), "UCERF3",
						SectionRecurrenceComparePlot.SectType.SUBSECTION, mapper, riMinMags);
				subSectCompare.initialize(getName(), outputDir, "interevent_sub_sects");
				plots.add(subSectCompare);
			}
			lines.add("### Element/Subsection Interevent Time Comparisons");
			lines.add("");
			lines.add("#### Element Interevent Time Comparisons");
			lines.add(topLink);
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.addLine("Min Mag", "Scatter", "2-D Hist");
			for (double riMinMag : riMinMags) {
				String prefix = "interevent_elements_m";
				if (riMinMag == Math.floor(riMinMag))
					prefix += (int)riMinMag;
				else
					prefix += (float)riMinMag;
				table.initNewLine();
				table.addColumn("**M≥"+(float)riMinMag+"**");
				table.addColumn("![Element Scatter]("+outputDir.getName()+"/"+prefix+"_scatter.png)");
				table.addColumn("![Element 2-D Hist]("+outputDir.getName()+"/"+prefix+"_hist2D.png)");
				table.finalizeLine();
			}
			lines.addAll(table.build());
			lines.add("");
			lines.add("#### Subsection Interevent Time Comparisons");
			lines.add(topLink);
			lines.add("");
			if (minFractForInclusion > 0) {
				lines.add("*Subsections participate in a rupture if at least "+(float)(minFractForInclusion*100d)+" % of its area ruptures*");
				lines.add("");
			}
			table = MarkdownUtils.tableBuilder();
			table.addLine("Min Mag", "Scatter", "2-D Hist");
			for (double riMinMag : riMinMags) {
				String prefix = "interevent_sub_sects_m";
				if (riMinMag == Math.floor(riMinMag))
					prefix += (int)riMinMag;
				else
					prefix += (float)riMinMag;
				table.initNewLine();
				table.addColumn("**M≥"+(float)riMinMag+"**");
				table.addColumn("![Subsection Scatter]("+outputDir.getName()+"/"+prefix+"_scatter.png)");
				table.addColumn("![Subsection 2-D Hist]("+outputDir.getName()+"/"+prefix+"_hist2D.png)");
				table.finalizeLine();
			}
			lines.addAll(table.build());
		}
		
		boolean ca = fm instanceof FaultModels || fm instanceof NSHM23_FaultModels;

		if (ca && plotsSet.contains(StandardPlots.PALEO_RECURRENCE) && hasFM) {
			File paleoCSVFile = new File(outputDir, "paleo_recurrence.csv");
			if (replot || !paleoCSVFile.exists()) {
				PaleoRecurrencePlot plot = new PaleoRecurrencePlot(getElements(), getSubSectMapper());
				plot.initialize(getName(), outputDir, "paleo_recurrence");
				plots.add(plot);
			}
			lines.add("");
			lines.add("### Paleo Recurrence Plots");
			lines.add(topLink); lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.addLine("![Paleo Plot]("+outputDir.getName()+"/paleo_recurrence_raw_sect_rate.png)",
					"![Paleo Plot]("+outputDir.getName()+"/paleo_recurrence_paleo_sect_rate.png)");
			table.addLine("![Paleo Plot]("+outputDir.getName()+"/paleo_recurrence_raw_elem_rate.png)",
					"![Paleo Plot]("+outputDir.getName()+"/paleo_recurrence_paleo_elem_rate.png)");
			lines.addAll(table.build());
			lines.add("");
			if (paleoCSVFile.exists()) {
				CSVFile<String> csv = CSVFile.readFile(new File(outputDir, "paleo_recurrence.csv"), true);
				lines.addAll(MarkdownUtils.tableFromCSV(csv, false).build());
			}
		}

		lines.add("");
		if (ca && plotsSet.contains(StandardPlots.PALEO_OPEN_INTERVAL)) {
			lines.add("### Paleo Open Interval Plots");
			lines.add(topLink); lines.add("");
			lines.add("#### Paleo Open Interval Plots, Biasi and Sharer 2019");
			lines.add(topLink); lines.add("");
			lines.add("These plots use the 5 paleoseismic sites identified in Biasi & Scharer (2019) on the Hayward, N. SAF, S. SAF, and SJC faults. "
					+ "By default, a rupture is counted at a paleo site if the nearest element (at the surface) slips any amount. We also alternatively "
					+ "apply a probability of detection model. Those results are marked as 'Prob. Filtered'.");
			lines.add("");
			lines.addAll(getPaleoPlotLines(outputDir, "paleo_open_biasi"));
			
			if (replot || !new File(outputDir, "paleo_open_biasi_prob.png").exists()) {
				PaleoOpenIntervalPlot paleoPlot = new PaleoOpenIntervalPlot(getElements(), PaleoOpenIntervalPlot.getSetBiasi2019(), 100);
				paleoPlot.initialize(getName(), outputDir, "paleo_open_biasi");
				plots.add(paleoPlot);
			}
			lines.add("");
			lines.add("#### Paleo Open Interval Plots, UCERF3");
			lines.add(topLink); lines.add("");
			lines.add("These plots use the full set of UCERF3 paleoseismic sites. "
					+ "By default, a rupture is counted at a paleo site if the nearest element (at the surface) slips any amount. We also alternativesly"
					+ "apply a probability of detection model. Those results are marked as 'Prob. Filtered'.");
			lines.add("");
			lines.addAll(getPaleoPlotLines(outputDir, "paleo_open_ucerf3"));
			
			if (replot || !new File(outputDir, "paleo_open_ucerf3_prob.png").exists()) {
				PaleoOpenIntervalPlot paleoPlot = new PaleoOpenIntervalPlot(getElements(), PaleoOpenIntervalPlot.getSetUCERF3(), 100);
				paleoPlot.initialize(getName(), outputDir, "paleo_open_ucerf3");
				plots.add(paleoPlot);
			}
		}
		
		if (plotsSet.contains(StandardPlots.MOMENT_RATE_VARIABILITY)) {
			if (getDurationYears() > MomentRateVaribilityPlot.MIN_CAT_LEN_YEARS) {
				if (replot || !new File(outputDir, "moment_variability_time_series.png").exists()) {
					MomentRateVaribilityPlot momVarPlot = new MomentRateVaribilityPlot(25, 2000);
					momVarPlot.initialize(getName(), outputDir, "moment_variability");
					plots.add(momVarPlot);
				}
				lines.add("");
				lines.add("### Moment Release Variability Plots");
				lines.add(topLink); lines.add("");
				lines.add("We first create a tapered moment release time series for the entire catalog. Each event's moment is distributed "
						+"across a 25 year Hanning (cosine) taper. Here is a plot of a random 2,000 year section of this time series:");
				lines.add("");
				lines.add("![Time Series]("+outputDir.getName()+"/moment_variability_time_series.png)");
				lines.add("");
				lines.add("We then compute Welch's power spectral density estimate on the entire time series. Results are plotted below, "
						+ "with a Poisson randomization of the catalog also plotted in a gray line, and the 95% confidence bounds from "
						+ MomentRateVaribilityPlot.num_poisson+" realizations as a light gray shaded area. Significant deviations outside the "
						+ "Poisson confidence intervals indicate synchronous behaviour.");
				lines.add("");
				lines.add("![Welch PSD]("+outputDir.getName()+"/moment_variability_welch.png)");
			}
		}
		
		if (plotsSet.contains(StandardPlots.ELASTIC_REBOUND_TRIGGERING) && hasFM) {
			double[] maxTimes = { 1d, 10d, 100d }; 
			double[] triggerMinMags = { 6d, 6.5d, 7d };
			if (replot || !new File(outputDir, "elastic_rebound_triggering_m7_1yr.png").exists()) {
				for (double triggerMinMag : triggerMinMags) {
					if (triggerMinMag < minMag)
						continue;
					ElasticReboundTriggeringPlot plot = new ElasticReboundTriggeringPlot(getSubSectMapper(), triggerMinMag, maxTimes);
					plot.initialize(getName(), outputDir, "elastic_rebound_triggering_m"+optionalDigitDF.format(triggerMinMag));
					plots.add(plot);
				}
			}
			lines.add("");
			lines.add("### Trigger Hypocenter Statistics Within Previous Rupture Area");
			lines.add(topLink); lines.add("");
			
			String[][] exampleArray = null;
			for (int i=0; i<triggerMinMags.length; i++) {
				double triggerMinMag = triggerMinMags[i];
				String prefix = "elastic_rebound_triggering_m"+optionalDigitDF.format(triggerMinMag);
				File example1 = new File(outputDir, prefix+"_example_re_rup.png");
				File example2 = new File(outputDir, prefix+"_example_new_hypo.png");
				if (example1.exists() || example2.exists()) {
					if (exampleArray == null)
						exampleArray = new String[2][triggerMinMags.length];
					exampleArray[0][i] = example1.exists() ? "![example]("+outputDir.getName()+"/"+example1.getName()+")" : "";
					exampleArray[1][i] = example2.exists() ? "![example]("+outputDir.getName()+"/"+example2.getName()+")" : "";
				}
			}
			if (exampleArray != null) {
				table = MarkdownUtils.tableBuilder();
				table.initNewLine();
				for (double triggerMinMag : triggerMinMags)
					if (triggerMinMag >= minMag)
						table.addColumn("M≥"+optionalDigitDF.format(triggerMinMag));
				table.finalizeLine();
				for (String[] vals : exampleArray)
					table.addLine(vals);
				lines.add("Example rupture plots:");
				lines.add("");
				lines.addAll(table.build());
				lines.add("");
			}
			
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			for (double triggerMinMag : triggerMinMags)
				if (triggerMinMag >= minMag)
					table.addColumn("M≥"+optionalDigitDF.format(triggerMinMag));
			table.finalizeLine();
				
			for (double maxTime : maxTimes) {
				table.initNewLine();
				for (double triggerMinMag : triggerMinMags) {
					if (triggerMinMag >= minMag) {
						String prefix = "elastic_rebound_triggering_m"+optionalDigitDF.format(triggerMinMag)
							+"_"+optionalDigitDF.format(maxTime)+"yr";
						table.addColumn("![hypocenter plot]("+outputDir.getName()+"/"+prefix+".png)");
					}
				}
				table.finalizeLine();
			}
			lines.addAll(table.build());
		}
		
		return plots;
	}
	
	private static String imageIfExists(File outputDir, String fileName, String title, String missingText) {
		if (new File(outputDir, fileName).exists())
			return "!["+title+"]("+outputDir.getName()+"/"+fileName+")";
		return missingText;
	}
	
	public List<String> writeStandardDiagnosticPlots(File outputDir, int skipYears, double minMag, boolean replot, String topLink,
			StandardPlots... standardPlots) throws IOException {
		List<String> lines = new ArrayList<>();
		HashSet<StandardPlots> plotsSet;
		if (standardPlots == null || standardPlots.length == 0) {
			plotsSet = new HashSet<>(EnumSet.allOf(StandardPlots.class));
		} else {
			plotsSet = new HashSet<>();
			for (StandardPlots plot : standardPlots)
				plotsSet.add(plot);
		}
		List<AbstractPlot> plots = buildStandardPlotLines(outputDir, skipYears, minMag, replot, topLink, lines, plotsSet);
		
		if (plots.isEmpty())
			return lines;
		
		for (AbstractPlot p : plots)
			p.setPlotSize(650, 600);
		
		Loader l = loader().minMag(minMag).skipYears(skipYears);
		
		System.out.println("Iterating through catalog and generating "+plots.size()+" plots");
		Iterable<RSQSimEvent> iterable = l.iterable();
		
		for (RSQSimEvent e : iterable)
			for (AbstractPlot p : plots)
				p.processEvent(e);
		
		System.out.println("Finalizing plots");
		for (AbstractPlot p : plots) {
//			System.out.println("Finalizing "+p.getClass());
			p.finalizePlot();
		}
		System.out.println("Done with plots");
		
		// lines could have changed, regenerate
		lines.clear();
		buildStandardPlotLines(outputDir, skipYears, minMag, false, topLink, lines, plotsSet);
		
		return lines;
	}
	
	private List<String> getPaleoPlotLines(File outputDir, String prefix) throws IOException {
		List<String> lines = new ArrayList<>();
		File sitesCSVFile = new File(outputDir, prefix+"_sites.csv");
		if (sitesCSVFile.exists()) {
			CSVFile<String> csv = CSVFile.readFile(sitesCSVFile, true);
			lines.add("**Paleoseismic sites table:**");
			lines.add("");
			lines.addAll(MarkdownUtils.tableFromCSV(csv, true).build());
			lines.add("");
		}
		lines.add("**Paleoseismic Plots:**");
		lines.add("");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		table.addColumn("![Count]("+outputDir.getName()+"/"+prefix+"_count.png)");
		table.addColumn("![Prob]("+outputDir.getName()+"/"+prefix+"_prob.png)");
		table.finalizeLine();
		lines.addAll(table.build());
		File probsCSVFile = new File(outputDir, prefix+"_probs.csv");
		if (sitesCSVFile.exists()) {
			CSVFile<String> csv = CSVFile.readFile(probsCSVFile, true);
			lines.add("");
			lines.add("**Open interval probabilities table:**");
			lines.add("");
			lines.addAll(MarkdownUtils.tableFromCSV(csv, true).build());
		}
		return lines;
	}

	@Override
	public Element toXMLMetadata(Element root) {
		Element el = root.addElement(XML_METADATA_NAME);
		
		el.addAttribute("name", name);
		el.addAttribute("author", author);
		el.addAttribute("dateMillis", date == null ? "0" : date.getTimeInMillis()+"");
		el.addAttribute("metadata", metadata);
		if (fm != null && fm instanceof Enum<?>) {
			el.addAttribute("fm", ((Enum<?>)fm).name());
			el.addAttribute("fmClass", fm.getClass().getName());
		}
		if (dm != null && dm instanceof Enum<?>) {
			el.addAttribute("dm", ((Enum<?>)dm).name());
			el.addAttribute("dmClass", dm.getClass().getName());
		}
		el.addAttribute("utmZone", utmZone+"");
		el.addAttribute("utmBand", utmBand+"");
		try {
			Map<Integer, Double> slipVels = getSlipVelocities();
			el.addAttribute("slipVel", constSlipVel+"");
			if (!Double.isFinite(constSlipVel)) {
				// write individual
				Element slipVelsEl = el.addElement("SlipVelocities");
				for (int patchID : slipVels.keySet()) {
					Element patchEl = slipVelsEl.addElement("Patch");
					patchEl.addAttribute("id", patchID+"");
					patchEl.addAttribute("velocity", slipVels.get(patchID)+"");
				}
			}
			el.addAttribute("aveArea", getAveArea()+"");
			el.addAttribute("numEvents", getNumEvents()+"");
			el.addAttribute("durationYears", getDurationYears()+"");
		} catch (Exception e) {}
		
		return root;
	}
	
	static RSQSimCatalog fromXMLMetadata(Element el) {
		String name = el.attributeValue("name");
		String author = el.attributeValue("author");
		long dateMillis = Long.parseLong(el.attributeValue("dateMillis"));
		GregorianCalendar cal = new GregorianCalendar();
		cal.setTimeInMillis(dateMillis);
		String metadata = el.attributeValue("metadata");
		RupSetFaultModel fm = null;
		if (el.attribute("fm") != null) {
			Attribute classAtt = el.attribute("fmClass");
			String enumName = el.attributeValue("fm");
			if (classAtt == null) {
				// assume U3
				fm = FaultModels.valueOf(enumName);
			} else {
				try {
					Class<? extends Enum<?>> fmClass = (Class<? extends Enum<?>>) Class.forName(classAtt.getValue());
					if (!fmClass.isEnum())
						fmClass = (Class<? extends Enum<?>>) fmClass.getEnclosingClass();
					for (Enum<?> eConst : fmClass.getEnumConstants()) {
						if (eConst.name().equals(enumName)) {
							fm = (RupSetFaultModel)eConst;
							break;
						}
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		RupSetDeformationModel dm = null;
		if (el.attribute("dm") != null) {
			Attribute classAtt = el.attribute("dmClass");
			String enumName = el.attributeValue("dm");
			if (classAtt == null) {
				// assume U3
				dm = DeformationModels.valueOf(enumName);
			} else {
				try {
					Class<? extends Enum<?>> dmClass = (Class<? extends Enum<?>>) Class.forName(classAtt.getValue());
					if (!dmClass.isEnum())
						dmClass = (Class<? extends Enum<?>>) dmClass.getEnclosingClass();
					for (Enum<?> eConst : dmClass.getEnumConstants()) {
						if (eConst.name().equals(enumName)) {
							dm = (RupSetDeformationModel)eConst;
							break;
						}
					}
				} catch (ClassNotFoundException | ClassCastException e) {
					e.printStackTrace();
				}
			}
		}
		double constSlipVel = Double.parseDouble(el.attributeValue("slipVel"));
		Map<Integer, Double> slipVels = null;
		if (!Double.isFinite(constSlipVel)) {
			Element slipVelEl = el.element("SlipVelocities");
			if (slipVelEl != null) {
				slipVels = new HashMap<>();
				for (Element patchEl : XMLUtils.getSubElementsList(slipVelEl)) {
					int patchID = Integer.parseInt(patchEl.attributeValue("id"));
					double patchVel = Double.parseDouble(patchEl.attributeValue("velocity"));
					slipVels.put(patchID, patchVel);
				}
			}
		}
		double aveArea = Double.NaN;
		if (el.attribute("aveArea") != null)
			aveArea = Double.parseDouble(el.attributeValue("aveArea"));
		int numEvents = -1;
		if (el.attribute("numEvents") != null)
			numEvents = Integer.parseInt(el.attributeValue("numEvents"));
		double durationYears = Double.NaN;
		if (el.attribute("durationYears") != null)
			durationYears = Double.parseDouble(el.attributeValue("durationYears"));
		
		int utmZone = 11;
		char utmBand = 'N';
		
		if (el.attribute("utmZone") != null)
			utmZone = Integer.parseInt(el.attributeValue("utmZone"));
		if (el.attribute("utmBand") != null)
			utmBand = el.attributeValue("utmBand").charAt(0);
		
		RSQSimCatalog cat = new RSQSimCatalog(name, author, cal, metadata, fm, dm, utmZone, utmBand);
		cat.aveArea = aveArea;
		cat.numEvents = numEvents;
		cat.durationYears = durationYears;
		cat.constSlipVel = constSlipVel;
		cat.slipVels = slipVels;
		return cat;
	}
	
	public static void writeCatalogsIndex(File dir) throws IOException, DocumentException {
		writeCatalogsIndex(dir, false, null, null);
	}
	
	public static void writeCatalogsIndex(File dir, boolean plotMulti, String baselineModelDirName,
			Map<String, List<String>> variationGroupings) throws IOException, DocumentException {
		// sort by date, newest first
		List<Long> times = new ArrayList<>();
		List<RSQSimCatalog> catalogs = new ArrayList<>();
		for (File subDir : dir.listFiles()) {
			if (!subDir.isDirectory())
				continue;
			File xmlFile = new File(subDir, "catalog.xml");
			if (!xmlFile.exists())
				continue;
			System.out.println("Writing catalog XML to "+xmlFile.getAbsolutePath());
			Document doc = XMLUtils.loadDocument(xmlFile);
			Element root = doc.getRootElement();
			Element el = root.element(XML_METADATA_NAME);
			RSQSimCatalog cat = fromXMLMetadata(el);
			times.add(cat.getDate().getTimeInMillis());
			cat.dir = subDir;
			catalogs.add(cat);
		}
		catalogs = ComparablePairing.getSortedData(times, catalogs);
		Collections.reverse(catalogs);
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Date", "Name", "Duration", "Element Area", "Description");
		for (RSQSimCatalog cat : catalogs) {
			table.initNewLine();
			
			table.addColumn(dateFormat.format(cat.getDate().getTime()));
			table.addColumn("["+cat.getName()+"]("+cat.dir.getName()+"#"+MarkdownUtils.getAnchorName(cat.getName())+")");
			try {
				table.addColumn(groupedIntDF.format(cat.getDurationYears())+" yrs");
			} catch (IOException e) {
				table.addColumn("");
			}
			try {
				table.addColumn(areaDF.format(cat.getAveArea())+" km");
			} catch (IOException e) {
				table.addColumn("");
			}
			table.addColumn(cat.getMetadata());
			
			table.finalizeLine();
		}
		
		List<String> lines = new LinkedList<>();
		lines.add("# RSQSim Catalogs Analysis");
		lines.add("");
		lines.addAll(table.build());
		
		if (plotMulti && catalogs.size() > 1)
			lines.addAll(writeMultiCatalogPlots(dir, catalogs, baselineModelDirName, variationGroupings));
		
		MarkdownUtils.writeReadmeAndHTML(lines, dir);
	}
	
	private static List<String> writeMultiCatalogPlots(File dir, List<RSQSimCatalog> catalogs, String baselineModelDirName,
			Map<String, List<String>> variationGroupings) throws IOException {
		Map<String, RSQSimCatalog> nameToCatMap = new HashMap<>();
		for (RSQSimCatalog catalog : catalogs)
			nameToCatMap.put(catalog.getName(), catalog);
		
		RSQSimCatalog baselineCatalog = null;
		if (baselineModelDirName != null)
			baselineCatalog = nameToCatMap.get(baselineModelDirName);
		
		Map<String, List<RSQSimCatalog>> variationCatalogGroupoings = null;
		if (variationGroupings != null) {
			variationCatalogGroupoings = new HashMap<>();
			for (String cat : variationGroupings.keySet()) {
				List<RSQSimCatalog> catalogsGroup = new ArrayList<>();
				for (String dirName : variationGroupings.get(cat))
					catalogsGroup.add(nameToCatMap.get(dirName));
				variationCatalogGroupoings.put(cat, catalogsGroup);
			}
		}
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(20);
		plotPrefs.setAxisLabelFontSize(22);
		plotPrefs.setPlotLabelFontSize(24);
		plotPrefs.setLegendFontSize(22);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		File resourcesDir = new File(dir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		System.out.println("Writing multi-catalog plots...");
		ArrayList<String> lines = new ArrayList<>();
		lines.add("## Multi-Catalog Plots");
		lines.add("");
		int tocIndex = lines.size();
		lines.add("");
		String topLink = "*[(top)](#"+MarkdownUtils.getAnchorName("Multi-Catalog Plots")+")*";
		if (baselineCatalog != null) {
			lines.add("Baseline catalog: ["+baselineCatalog.getName()+"]("+baselineCatalog.dir.getName()
				+"#"+MarkdownUtils.getAnchorName(baselineCatalog.getName())+")");
			lines.add("");
		}
		
		// MFD plots
		Map<RSQSimCatalog, DiscretizedFunc> mfdMap = new HashMap<>();
		
		for (RSQSimCatalog catalog : catalogs) {
			File catalogDir = catalog.getCatalogDir();
			File csvFile = new File(catalogDir, "resources/mfd.csv");
			if (csvFile.exists()) {
				CSVFile<String> csv = CSVFile.readFile(csvFile, true);
				DiscretizedFunc mfd = new ArbitrarilyDiscretizedFunc();
				for (int row=1; row<csv.getNumRows(); row++) {
					double mag = Double.parseDouble(csv.get(row, 0));
					double cumRate = Double.parseDouble(csv.get(row, 2));
					mfd.set(mag, cumRate);
				}
				mfdMap.put(catalog, mfd);
			}
		}
		
		if (mfdMap.size() > 1) {
			boolean hasBaseline = baselineCatalog != null && mfdMap.get(baselineCatalog) != null;
			
			DiscretizedFunc baselineFunc = null;
			if (hasBaseline) {
				baselineFunc = mfdMap.get(baselineCatalog);
				baselineFunc.setName(baselineCatalog.getName());
			}
			
			MFDPlot.plotMultiMFDs(mfdMap.values(), baselineFunc, resourcesDir, "mfds");
			
			lines.add("### MFDs");
			lines.add(topLink); lines.add("");
			lines.add("![MFDs]("+resourcesDir.getName()+"/mfds.png)");
		}
		
		if (variationGroupings != null) {
			if (mfdMap.isEmpty()) {
				lines.add("### MFDs");
				lines.add(topLink); lines.add("");
			} else {
				lines.add("");
				lines.add("#### MFD Variations Table");
				lines.add("");
			}
			
			lines.addAll(buildVariationTable(dir, variationCatalogGroupoings, "mfd.png", 4));
		}
		
		// Mag-Area
		Map<RSQSimCatalog, CSVFile<String>> maMap = new HashMap<>();

		for (RSQSimCatalog catalog : catalogs) {
			File catalogDir = catalog.getCatalogDir();
			File csvFile = new File(catalogDir, "resources/mag_area.csv");
			if (csvFile.exists()) {
				CSVFile<String> csv = CSVFile.readFile(csvFile, true);
				maMap.put(catalog, csv);
			}
		}

		if (maMap.size() > 1) {
			boolean hasBaseline = baselineCatalog != null && maMap.get(baselineCatalog) != null;

			CSVFile<String> baselineCSV = null;
			String baselineName = null;
			if (hasBaseline) {
				baselineCSV = maMap.get(baselineCatalog);
				baselineName = baselineCatalog.getName();
			}
			
			MagAreaScalingPlot.plotMultiMagArea(maMap.values(), baselineCSV, baselineName, false, true,
					new double[] {0.025, 0.975}, resourcesDir, "mag_areas");

			lines.add("### Magnitude-Area Plots");
			lines.add(topLink); lines.add("");
			lines.add("![Mag Areas]("+resourcesDir.getName()+"/mag_areas.png)");
		}
		
		if (variationGroupings != null) {
			if (maMap.isEmpty()) {
				lines.add("### Magnitude-Area Plots");
				lines.add(topLink); lines.add("");
			} else {
				lines.add("");
				lines.add("#### M-A Variations Table");
				lines.add("");
			}
			
			lines.addAll(buildVariationTable(dir, variationCatalogGroupoings, "mag_area_hist2D.png", 4));
		}
		
		// Rupture Velocity
		Map<RSQSimCatalog, CSVFile<String>> magVelMap = new HashMap<>();

		for (RSQSimCatalog catalog : catalogs) {
			File catalogDir = catalog.getCatalogDir();
			File csvFile = new File(catalogDir, "resources/rupture_velocity_scatter_fractiles.csv");
			if (csvFile.exists()) {
				CSVFile<String> csv = CSVFile.readFile(csvFile, true);
				magVelMap.put(catalog, csv);
			}
		}

		if (magVelMap.size() > 1) {
			boolean hasBaseline = baselineCatalog != null && magVelMap.get(baselineCatalog) != null;

			CSVFile<String> baselineCSV = null;
			String baselineName = null;
			if (hasBaseline) {
				baselineCSV = magVelMap.get(baselineCatalog);
				baselineName = baselineCatalog.getName();
			}

			RuptureVelocityPlot.plotMultiMagVels(magVelMap.values(), baselineCSV, baselineName, false, false,
					new double[] {0.025, 0.975}, resourcesDir, "mag_velocities");

			lines.add("### Rupture Velocity vs Magnitude");
			lines.add(topLink); lines.add("");
			lines.add("![Velocities]("+resourcesDir.getName()+"/mag_velocities.png)");
		}
		
		if (variationGroupings != null) {
			if (magVelMap.isEmpty()) {
				lines.add("### Rupture Velocity vs Magnitude");
				lines.add(topLink); lines.add("");
			} else {
				lines.add("");
				lines.add("#### Rupture Velocity vs Magnitude Variations Table");
				lines.add("");
			}
			
			lines.addAll(buildVariationTable(dir, variationCatalogGroupoings, "rupture_velocity_scatter.png", 4));
		}
		
		Map<RSQSimCatalog, CSVFile<String>> distVelMap = new HashMap<>();

		for (RSQSimCatalog catalog : catalogs) {
			File catalogDir = catalog.getCatalogDir();
			File csvFile = new File(catalogDir, "resources/rupture_velocity_vs_dist.csv");
			if (csvFile.exists()) {
				CSVFile<String> csv = CSVFile.readFile(csvFile, true);
				distVelMap.put(catalog, csv);
			}
		}

		if (distVelMap.size() > 1) {
			boolean hasBaseline = baselineCatalog != null && distVelMap.get(baselineCatalog) != null;

			CSVFile<String> baselineCSV = null;
			String baselineName = null;
			if (hasBaseline) {
				baselineCSV = distVelMap.get(baselineCatalog);
				baselineName = baselineCatalog.getName();
			}

			RuptureVelocityPlot.plotMultiDistVels(distVelMap.values(), baselineCSV, baselineName, resourcesDir, "dist_velocities");

			lines.add("### Rupture Velocity vs Distance");
			lines.add(topLink); lines.add("");
			lines.add("![Velocities]("+resourcesDir.getName()+"/dist_velocities.png)");
		}
		
		if (variationGroupings != null) {
			if (distVelMap.isEmpty()) {
				lines.add("### Rupture Velocity vs Distance");
				lines.add(topLink); lines.add("");
			} else {
				lines.add("");
				lines.add("#### Rupture Velocity vs Distance Variations Table");
				lines.add("");
			}
			
			lines.addAll(buildVariationTable(dir, variationCatalogGroupoings, "rupture_velocity_vs_dist.png", 4));
		}
		
		if (variationGroupings != null) {
			for (double intereventMags : new double[] { 7 }) {
				lines.add("### M"+optionalDigitDF.format(intereventMags)+" Element Interevent Time Comparisons");
				lines.add(topLink); lines.add("");
				lines.addAll(buildVariationTable(dir, variationCatalogGroupoings,
						"interevent_elements_m"+optionalDigitDF.format(intereventMags)+"_hist2D.png", 4));
				lines.add("");
				lines.add("### M"+optionalDigitDF.format(intereventMags)+" Subsection Interevent Time Comparisons");
				lines.add(topLink); lines.add("");
				lines.addAll(buildVariationTable(dir, variationCatalogGroupoings,
						"interevent_sub_sects_m"+optionalDigitDF.format(intereventMags)+"_hist2D.png", 4));
			}
		}
		
		// Palo Probs
		if (variationGroupings != null) {
			lines.add("### Paleo Open Interval Plots");
			lines.add(topLink); lines.add("");
			lines.add("#### Paleo Open Interval Plots, Biasi and Sharer 2019");
			lines.add(topLink); lines.add("");
			lines.addAll(buildVariationTable(dir, variationCatalogGroupoings, "paleo_open_biasi_prob.png", 4));
			lines.add("");
			lines.add("#### Paleo Open Interval Plots, UCERF3");
			lines.add(topLink); lines.add("");
			lines.addAll(buildVariationTable(dir, variationCatalogGroupoings, "paleo_open_ucerf3_prob.png", 4));
		}
		
		// Welch PSD
		if (variationGroupings != null) {
			lines.add("### Moment Release Variability Welch PSDs");
			lines.add(topLink); lines.add("");
			lines.addAll(buildVariationTable(dir, variationCatalogGroupoings, "moment_variability_welch.png", 4));
		}
		
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		
		System.out.println("DONE");
		
		return lines;
	}
	
	private static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");
	
	private static List<String> buildVariationTable(File dir, Map<String, List<RSQSimCatalog>> variationGroupings,
			String fileName, int wrapCols) {
		List<String> lines = new ArrayList<>();
		for (String groupName : variationGroupings.keySet()) {
			TableBuilder table = MarkdownUtils.tableBuilder();
			List<String> images = new ArrayList<>();
			
			table.initNewLine();
			for (RSQSimCatalog catalog : variationGroupings.get(groupName)) {
				File catDir = catalog.getCatalogDir();
				File resourcesDir = new File(catDir, "resources");
				File imageFile = new File(resourcesDir, fileName);
				if (!imageFile.exists())
					continue;
				table.addColumn("<p align=\"center\">**"+catalog.getMetadata()+"**</p>");
				images.add(catDir.getName()+"/resources/"+fileName);
			}
			table.finalizeLine();
			
			if (images.isEmpty())
				continue;
			
			table.initNewLine();
			for (String image : images)
				table.addColumn("![plot]("+image+")");
			table.finalizeLine();
			
			table.wrap(wrapCols, 0);
			lines.add("**"+groupName+"**");
			lines.add("");
			lines.addAll(table.build());
			lines.add("");
		}
		
		return lines;
	}
	
	static class CatEnumDateComparator implements Comparator<Catalogs> {

		@Override
		public int compare(Catalogs o1, Catalogs o2) {
			// reverse sorted, newest first
			return o2.catalog.getDate().compareTo(o1.catalog.getDate());
		}
		
	}
	
	private static final File[] catalogLocations;
	static {
		catalogLocations = new File[] {
				// USC HPC
				new File("/project/scec_608/rsqsim/catalogs/kmilner"),
				new File("/project/scec_608/rsqsim/catalogs/shaw"),
				new File("/project/scec_608/rsqsim/catalogs/gilchrij"),
				new File("/project/scec_608/rsqsim/catalogs/gilchrij/cybershake"),
				new File("/project/scec_608/rsqsim/catalogs/gilchrij/paramSweep"),
				// Kevin's laptop
				new File("/data/kevin/simulators/catalogs"),
				new File("/data/kevin/simulators/catalogs/bruce"),
				// TACC Stampede2
				new File("/work/00950/kevinm/stampede2/simulators/catalogs")
		};
	}
	
	public static File locateCatalog(String dirName, String... requiredFiles) {
		catalogDir:
		for (File dir : catalogLocations) {
			File catalogDir = new File(dir, dirName);
			if (!catalogDir.exists())
				continue;
			for (String requriedFile : requiredFiles)
				if (!new File(catalogDir, requriedFile).exists())
					continue catalogDir;
			return catalogDir;
		}
		return null;
	}
	
	public static void main(String args[]) throws IOException, DocumentException {
		File gitDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");
		
		boolean overwriteIndividual = true;
		boolean replot = true;
		
//		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		Catalogs[] cats = Catalogs.values();
		Arrays.sort(cats, new CatEnumDateComparator());
		// new catalogs
//		GregorianCalendar minDate = cal(2021, 10, 1);
		GregorianCalendar minDate = cal(2024, 8, 1);
		for (Catalogs cat : cats) {
		// specific catalog
//		GregorianCalendar minDate = cal(2000, 1, 1);
//		for (Catalogs cat : new Catalogs[] {
//////				Catalogs.BRUCE_4983_STITCHED,
//////////				Catalogs.BRUCE_2585,
//////////				Catalogs.BRUCE_2585_1MYR,
//////////				Catalogs.BRUCE_2740,
//////////				Catalogs.BRUCE_3062,
//////////				Catalogs.BRUCE_4860,
//////////				Catalogs.JG_tunedBase1m_ddotEQmod,
//////////				Catalogs.JG_tuneBase1m,
//				Catalogs.BRUCE_5672
//				}) {
		// all catalogs
//		GregorianCalendar minDate = cal(2000, 1, 1);
//		for (Catalogs cat : cats) {
			
			if (cat.catalog.getDate().before(minDate))
				continue;
			RSQSimCatalog catalog = cat.instance();
			System.out.print(catalog.getName()+" ? ");
			File catGitDir = new File(gitDir, catalog.getCatalogDir().getName());
			Preconditions.checkState(catGitDir.exists() || catGitDir.mkdir());
			File xmlFile = new File(catGitDir, "catalog.xml");
			if (xmlFile.exists())
				System.out.println("exists");
			else
				System.out.println("missing");
			if (overwriteIndividual || !xmlFile.exists()) {
				System.out.println("\twriting summary");
				catalog.writeMarkdownSummary(catGitDir, true, replot);
			}
		}
		
		writeCatalogsIndex(gitDir);
	}
}
