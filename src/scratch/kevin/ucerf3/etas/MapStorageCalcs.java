package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.utils.FaultSystemIO;

public class MapStorageCalcs {

	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemRupSet rupSet = FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2016_06_15-haywired_m7-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-combined/results_m5_preserve.bin");
		List<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(catalogFile);
		
		int numGridded = 0;
		int numFSS = 0;
		
		double cutoff = 200d;
		double spacing = 0.1;
		
		long pointsGridded = 0l;
		long pointsFSS = 0l;
		
		HashSet<Integer> griddedIndexes = new HashSet<Integer>();
		HashSet<Integer> fssIndexes = new HashSet<Integer>();
		
		for (int i=0; i<catalogs.size(); i++) {
			if (i % 1000 == 0) System.out.println("Processing catalog "+i+"/"+catalogs.size());
			List<ETAS_EqkRupture> catalog = catalogs.get(i);
			for (ETAS_EqkRupture rup : catalog) {
				if (rup.getFSSIndex() >= 0) {
					// supra seismogenic
					numFSS++;
					int fssIndex = rup.getFSSIndex();
					if (!fssIndexes.contains(fssIndex)) {
						// new rupture
						fssIndexes.add(fssIndex);
//						LocationList trace = rupSet.getSurfaceForRupupture(fssIndex, 1d, false).getUpperEdge();
//						GriddedRegion reg = new GriddedRegion(trace, cutoff, spacing, null);
						
//						List<LocationList> tracesByParent = Lists.newArrayList();
//						int parentID = -1;
//						List<String> parentNames = Lists.newArrayList();
//						List<List<String>> sectNames = Lists.newArrayList();
//						for (int s : rupSet.getSectionsIndicesForRup(fssIndex)) {
//							FaultSectionPrefData fsd = rupSet.getFaultSectionData(s);
//							int myParent = fsd.getParentSectionId();
//							if (myParent != parentID) {
//								tracesByParent.add(new LocationList());
//								parentID = myParent;
//								parentNames.add(fsd.getParentSectionName());
//								sectNames.add(new ArrayList<String>());
//							}
//							LocationList curList = tracesByParent.get(tracesByParent.size()-1);
//							for (Location loc : fsd.getFaultTrace()) {
//								if (curList.isEmpty() || !curList.last().equals(loc))
//									curList.add(loc);
//							}
////							tracesByParent.get(tracesByParent.size()-1).addAll(fsd.getFaultTrace());
//							sectNames.get(sectNames.size()-1).add(fsd.getSectionName());
//						}
//						List<Region> parentRegionsComposite = Lists.newArrayList();
//						for (int j=0; j<tracesByParent.size(); j++) {
//							LocationList parentTrace = tracesByParent.get(j);
////							List<Region> subRegions = regionsByParent.get(j);
////							Region compositeRegion = null;
////							for (int k=0; k<subRegions.size(); k++) {
////								Region subRegion = subRegions.get(k);
////								if (compositeRegion == null)
////									compositeRegion = new Region(subRegion);
////								else
////									compositeRegion = Region.union(compositeRegion, subRegion);
////								if (compositeRegion == null) {
////									System.out.println("Rupture has "+regionsByParent.size()+" parents");
////									System.out.println("\tFailed on: "+sectNames.get(j).get(k));
////									for (String parentName : parentNames)
////										System.out.println("\t"+Joiner.on("\n\t").join(sectNames.get(j)));
////								}
////								Preconditions.checkNotNull(compositeRegion);
////							}
//							System.out.println("Processing "+parentNames.get(j));
//							Region compositeRegion = new Region(parentTrace, cutoff);
//							parentRegionsComposite.add(compositeRegion);
//						}
						
//						CompoundSurface surf = (CompoundSurface) rupSet.getSurfaceForRupupture(fssIndex, 1d, false);
						List<List<FaultSectionPrefData>> sectsByParent = Lists.newArrayList();
						List<String> parents = Lists.newArrayList();
						int curParentID = -1;
						for (int s : rupSet.getSectionsIndicesForRup(fssIndex)) {
							FaultSectionPrefData fsd = rupSet.getFaultSectionData(s);
							if (fsd.getParentSectionId() != curParentID) {
								curParentID = fsd.getParentSectionId();
								sectsByParent.add(new ArrayList<FaultSectionPrefData>());
								parents.add(fsd.getParentSectionName());
							}
							sectsByParent.get(sectsByParent.size()-1).add(fsd);
						}
						
						List<Region> parentRegionsComposite = Lists.newArrayList();
						for (List<FaultSectionPrefData> sects : sectsByParent) {
							// build composite surface for just this parent section
							List<RuptureSurface> surfs = Lists.newArrayList();
							for(FaultSectionPrefData fltData : sects)
								surfs.add(fltData.getStirlingGriddedSurface(spacing, false, true));
							RuptureSurface compound;
							if (surfs.size() == 1)
								compound = surfs.get(0);
							else
								compound = new CompoundSurface(surfs);
							LocationList trace = compound.getEvenlyDiscritizedUpperEdge();
							// this can have duplicates in it, remove those
							for (int p=trace.size(); --p>0;) {
								Location p1 = trace.get(p);
								Location p2 = trace.get(p-1);
								double dist = LocationUtils.horzDistanceFast(p1, p2);
								if (dist < spacing)
									trace.remove(p);
//								if ((float)p1.getLatitude() == (float)p2.getLatitude()
//										&& (float)p1.getLongitude() == (float)p2.getLongitude())
//									trace.remove(p)
							}
							parentRegionsComposite.add(new Region(trace, cutoff));
						}
						
						Region compositeRegion = null;
						for (int j=0; j<parentRegionsComposite.size(); j++) {
							Region subRegion = parentRegionsComposite.get(j);
							Region prevRegion = compositeRegion;
							if (compositeRegion == null)
								compositeRegion = new Region(subRegion);
							else
								compositeRegion = Region.union(compositeRegion, subRegion);
							if (compositeRegion == null) {
								System.out.println("WTF?");
								double minDist = Double.POSITIVE_INFINITY;
								for (Location p1 : prevRegion.getBorder())
									for (Location p2 : subRegion.getBorder())
										minDist = Math.min(minDist, LocationUtils.horzDistanceFast(p1, p2));
								System.out.println("Min distance between prev region and new region: "+minDist);
								for (int k=0; k<parentRegionsComposite.size(); k++) {
									System.out.print("\t"+parents.get(k));
									if (k == j)
										System.out.println(" ********");
									else
										System.out.println();
								}
							}
							Preconditions.checkNotNull(compositeRegion);
						}
						Preconditions.checkNotNull(compositeRegion);
						GriddedRegion reg = new GriddedRegion(compositeRegion, spacing, null);
						pointsFSS += reg.getNodeCount();
					}
				} else {
					numGridded++;
					if (!griddedIndexes.contains(rup.getNthERF_Index())) {
						// new rupture
						griddedIndexes.add(rup.getNthERF_Index());
						GriddedRegion reg = new GriddedRegion(rup.getHypocenterLocation(), cutoff, spacing, null);
						pointsGridded += reg.getNodeCount();
					}
				}
			}
		}
		System.out.println(fssIndexes.size()+" unique FSS ruptures ("+numFSS+" total)");
		System.out.println(pointsFSS+" FSS grid nodes to be stored");
		printDiskUsage(pointsFSS);
		
		System.out.println(griddedIndexes.size()+" unique gridded ruptures ("+numGridded+" total)");
		System.out.println(pointsGridded+" FSS grid nodes to be stored");
		printDiskUsage(pointsGridded);
	}
	
	private static void printDiskUsage(long numPoints) {
		long valsToSave = numPoints*2l;
		long bytesFloat = valsToSave*4l;
		long bytesDouble = valsToSave*8l;
		double gbFloat = (double)bytesFloat/1073741824d;
		double gbDouble = (double)bytesDouble/1073741824d;
		System.out.println("\t"+valsToSave+" unique values");
		System.out.println("\t"+(float)gbFloat+" GB as floats");
		System.out.println("\t"+(float)gbDouble+" GB as doubles");
	}

}
