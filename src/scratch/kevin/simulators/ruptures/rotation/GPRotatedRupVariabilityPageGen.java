package scratch.kevin.simulators.ruptures.rotation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.MathArrays;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_Wrappers;
import org.opensha.sha.imr.mod.impl.BaylessSomerville2013DirectivityModifier;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.SRF_PointData;
import org.opensha.sha.simulators.utils.RSQSimSubSectEqkRupture;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SlipAlongSectAlgorithm;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.ruptures.ASK_EventData;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.RotationSpec;

public abstract class GPRotatedRupVariabilityPageGen extends RotatedRupVariabilityPageGen<GPRotatedRupture> {

	private List<GPRotatedRupture> ruptures;
	
	public GPRotatedRupVariabilityPageGen(GPRotatedRupVariabilityConfig config,
			double mag, SimulationRotDProvider<RotationSpec> prov, double[] calcPeriods) {
		super(config, mag, prov, calcPeriods);
		ruptures = config.getRuptureList();
	}
	
	@Override
	protected String getModelName() {
		return "Graves & Pitarka (2016)";
	}
	
	private Location calcHypocenter(List<SRF_PointData> srf) {
		double minTime = Double.POSITIVE_INFINITY;
		Location hypo = null;
		for (SRF_PointData point : srf) {
			if (point.getStartTime() < minTime) {
				minTime = point.getStartTime();
				hypo = point.getLocation();
			}
		}
		return hypo;
	}
	
	@Override
	protected EqkRupture buildGMPE_Rupture(GPRotatedRupture event) {
		double mag = event.src.getMag();
		double aveRake = event.src.getFocalMechanism().getRake();
		RuptureSurface surf = event.src.getSurface().getQuadSurface();
		Location hypo = calcHypocenter(event.srf);
		return new EqkRupture(mag, aveRake, surf, hypo);
	}
	
	@Override
	protected synchronized GPRotatedRupture getEvent(int eventID) {
		return ruptures.get(eventID);
	}
	
	@Override
	protected double calcVprop(GPRotatedRupture event) {
		double minTime = Double.POSITIVE_INFINITY;
		Location hypo = null;
		for (SRF_PointData rec : event.srf) {
			if (rec.getStartTime() < minTime) {
				minTime = rec.getStartTime();
				hypo = rec.getLocation();
			}
		}

		List<Double> vels = new ArrayList<>();
		for (SRF_PointData rec : event.srf) {
			if (rec.getLocation() == hypo)
				continue;
			double dist = LocationUtils.linearDistanceFast(hypo, rec.getLocation());
			
			double tDelta = rec.getStartTime() - minTime;
			if (tDelta == 0)
				continue;
			double vel = dist/(tDelta);
			vels.add(vel);
		}
		
		return DataUtils.median(Doubles.toArray(vels));
	}
	
	protected double getMag(GPRotatedRupture event) {
		return event.src.getMag();
	}
	
	protected double getArea(GPRotatedRupture event) {
		BBP_PlanarSurface surf = event.src.getSurface();
		double area = surf.getLength() * surf.getWidth();
		return area * 1e6; // km^2 to m^2
	}
	
	protected double getMaxSlip(GPRotatedRupture event) {
		double max = 0d;
		for (SRF_PointData point : event.srf)
			max = Math.max(max, point.getTotalSlip());
		return max;
	}
	
	protected double getMeanSlip(GPRotatedRupture event) {
		double mean = 0d;
		double totArea = 0d;
		for (int i=0; i<event.srf.size(); i++) {
			SRF_PointData point = event.srf.get(i);
			double area = point.getArea();
			mean += area*point.getTotalSlip();;
			totArea += area;
		}
		mean /= totArea;
		return mean;
	}
	
	protected double getSlipStdDev(GPRotatedRupture event) {
		double[] slips = new double[event.srf.size()];
		double[] weights = new double[slips.length];
		for (int i=0; i<slips.length; i++) {
			SRF_PointData point = event.srf.get(i);
			weights[i] = point.getArea();
			slips[i] = point.getTotalSlip();
		}
		double var = new Variance().evaluate(slips, MathArrays.normalizeArray(weights, weights.length));
		return Math.sqrt(var);
	}
	
	protected double getMeanMidSeisSlip(GPRotatedRupture event) {
		return Double.NaN;
	}
	
	@Override
	protected void plotExample(File resourcesDir, String prefix, double distance, List<Quantity> variedQuantities)
			throws IOException {
		List<Site> sites = new ArrayList<>();
		
		Double minMag = Double.POSITIVE_INFINITY;
		for (Double mag : magEventIDs.keySet())
			minMag = Double.min(minMag, mag);
		GPRotatedRupture exampleRupture = getEvent(magEventIDs.get(minMag).get(0));
		sites.add(this.sites.get(0));
		List<GPRotatedRupture> ruptures = new ArrayList<>();
		ruptures.add(exampleRupture);
		int numSourceAz = variedQuantities.contains(Quantity.SOURCE_AZIMUTH) ? numExampleRotations : 1;
		int numSiteToSourceAz = variedQuantities.contains(Quantity.SITE_TO_SOURTH_AZIMUTH) ? numExampleRotations : 1;
		GPRotatedRupVariabilityConfig config = new GPRotatedRupVariabilityConfig(sites, ruptures, new double[] {distance},
				numSourceAz, numSiteToSourceAz);
		
		config.plotRotations(resourcesDir, prefix, config.getRotations(), true);
	}
	
	@Override
	protected Table<Float, Double, XY_DataSet> calcDirectivityComparisons(Double mag, double[] periods,
			Site[] sites, File resourcesDir) {
		GPRotatedRupVariabilityConfig config = (GPRotatedRupVariabilityConfig)magConfigs.get(mag);
		
		List<Integer> eventIDs = config.getValues(Integer.class, Quantity.EVENT_ID);
		
		int[] periodIndexes = new int[periods.length];
		for (int i=0; i<periods.length; i++)
			periodIndexes[i] = Doubles.indexOf(calcPeriods, periods[i]);
		
		BaylessSomerville2013DirectivityModifier bs = new BaylessSomerville2013DirectivityModifier();
		
		Table<Float, Double, XY_DataSet> distPeriodXYs = HashBasedTable.create();
		
		TableBuilder table = null;
		File dirDir = null;
		if (DIRECTIVITY_DEBUG) {
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			for (double period : periods)
				table.addColumn(optionalDigitDF.format(period)+"s Map Plot");
			table.addColumn("Residuals Plot");
			table.finalizeLine();
			dirDir = new File(resourcesDir, DIRECTIVITY_DEBUG_DIRNAME);
		}
		
		try {
			for (int i=0; i<eventIDs.size(); i++) {
				// build simple rupture description
				int eventID = eventIDs.get(i);
				
				// rotated such that strike is zero, e.g. oriented facing north
				// aki & richards, dips to right (if dips)
				GPRotatedRupture event = config.getInitiallyOrientedRupture(eventID);
				
				Location centroid = event.centroid;
				Location hypo = calcHypocenter(event.srf);
				System.out.println("Hypocenter: "+hypo);
				RuptureSurface surf = event.src.getSurface().getQuadSurface();
				double aveRake = event.src.getFocalMechanism().getRake();
				EqkRupture rup = new EqkRupture(mag, aveRake, surf, hypo);

				Table<Float, Float, Location> azDistLocMap = HashBasedTable.create();
				Table<Float, Float, double[]> azDistFdMap = HashBasedTable.create();
				for (Float sourceAz : config.getValues(Float.class, Quantity.SOURCE_AZIMUTH)) {
					float siteAz = 180 - sourceAz;
					for (Float distance : distances) {
						Location siteLoc = LocationUtils.location(centroid, Math.toRadians(siteAz), distance);
						// correct distance
						for (int j=0; j<3; j++) {
							double minDist = BBP_PartBValidationConfig.DIST_JB ?
									surf.getDistanceJB(siteLoc) : surf.getDistanceRup(siteLoc);
							// positive means we're too close and need to move further
							double distDiff = distance - minDist;
							siteLoc = LocationUtils.location(siteLoc, Math.toRadians(siteAz), distDiff);
						}
						
						double[] fds = new double[periods.length];
						for (int p=0; p<periods.length; p++)
							fds[p] = bs.getFd(rup, siteLoc, periods[p]);
						
						azDistFdMap.put(sourceAz, distance, fds);
						azDistLocMap.put(sourceAz, distance, siteLoc);
					}
				}
				
				Map<RotationSpec, double[]> diffsMap = new HashMap<>();
				
				for (Float distance : distances) {
					EventTerm eventTerm = eventTermCache.get(
							new EventTermKey(eventID, mag, distance, sites));
					
					List<RotationSpec> rotations = eventTerm.rotations;
					for (int j=0; j<rotations.size(); j++) {
						RotationSpec rot = rotations.get(j);
						
						float sourceAz = nullAsZero(rot.sourceAz);
						
						double[] fds = azDistFdMap.get(sourceAz, distance);
						
						double[] diffs = new double[periods.length];
						
						for (int p=0; p<periods.length; p++) {
							double meanVal = eventTerm.eventTerms[periodIndexes[p]];
							double value = eventTerm.rotLogVals.get(j)[periodIndexes[p]];
							diffs[p] = value - meanVal;
							
							XY_DataSet xy = distPeriodXYs.get(distance, periods[p]);
							if (xy == null) {
								xy = new DefaultXY_DataSet();
								distPeriodXYs.put(distance, periods[p], xy);
							}
							xy.set(fds[p], diffs[p]);
						}
						diffsMap.put(rot, diffs);
					}
					
				}
				if (DIRECTIVITY_DEBUG) {
					File[] files = plotIndvRupDirectivity(eventID, null, rup, azDistFdMap, azDistLocMap,
							diffsMap, periods, dirDir);
					table.initNewLine();
					for (File file : files)
						table.addColumn("![plot]("+file.getName()+")");
					table.finalizeLine();
				}
			}
		} catch (Exception e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		if (DIRECTIVITY_DEBUG) {
			List<String> lines = new ArrayList<>();
			lines.add("# Directivity Debug");
			lines.add("");
			lines.add("Map plots show the re-oriented (with strike=0) raw rupture in black, "
					+ "the linear GMPE source used to predict directivity with a dashed line, "
					+ "and hypocenter with a star. Map background is empirical fD estimates, and mean "
					+ "simulated values at 20 and 50 km are overlaid. The chart on the right shows "
					+ "the simulated values for 20 km at each spectral period.");
			lines.add("");
			lines.addAll(table.build());
			lines.add("");
			File outputDir = new File(resourcesDir, DIRECTIVITY_DEBUG_DIRNAME);
			try {
				MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return distPeriodXYs;
	}

}
