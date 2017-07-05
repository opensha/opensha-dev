package scratch.kevin.simulators.momRateVariation;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.SimulatorElement;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.simulators.momRateVariation.UCERF3ComparisonAnalysis.UCERF3EventRecord;

public class UCERF3_ETASComparisons {
	
	public static List<List<SimulatorEvent>> loadUCERF3EtasCatalogs(List<List<ETAS_EqkRupture>> catalogs, FaultSystemSolution sol,
			Region region, Map<Integer, SimulatorElement> elems)
					throws IOException {
		
		List<List<SimulatorEvent>> eventsList = Lists.newArrayList();
		
		for (List<ETAS_EqkRupture> catalog : catalogs)
			eventsList.add(loadETASCatalogAsFakeSimEvents(sol, region, catalog, elems));
		
		return eventsList;
	}
	
	private static List<SimulatorEvent> loadETASCatalogAsFakeSimEvents(FaultSystemSolution sol, Region region,
			List<ETAS_EqkRupture> catalog, Map<Integer, SimulatorElement> elems) throws IOException {
		List<SimulatorEvent> events = Lists.newArrayList();
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		HashSet<Integer> sectIndexesInRegion;
		if (region == null) {
			sectIndexesInRegion = null;
		} else {
			sectIndexesInRegion = new HashSet<Integer>();
			for (FaultSectionPrefData sect : rupSet.getFaultSectionDataList()) {
				for (Location loc : sect.getFaultTrace()) {
					if (region.contains(loc)) {
						sectIndexesInRegion.add(sect.getSectionId());
						break;
					}
				}
			}
		}
		
		long startTime = catalog.get(0).getOriginTime();
		
		for (ETAS_EqkRupture rup : catalog) {
			if (rup.getFSSIndex() < 0)
				continue;
			long millis = rup.getOriginTime() - startTime;
			double secs = (double)millis / 1000d;
			
			EventRecord rec = new UCERF3EventRecord(elems, rupSet, rup.getFSSIndex(), secs);
			
			SimulatorEvent e = new SimulatorEvent(rec);
			
			events.add(e);
		}
		
		return events;
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File catalogsFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2015_11_30-spontaneous-1000yr-FelzerParams-mc20-full_td-noApplyLTR/results_m4.bin");
		File fssFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		
		File outputDir = new File("/home/kevin/Simulators/time_series/ucerf3_etas");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(catalogsFile);
		
		Region region = new CaliforniaRegions.RELM_SOCAL();
		
		FaultSystemSolution sol = FaultSystemIO.loadSol(fssFile);
		
		int[] windowLens = { 10, 25, 50, 75, 100, 150, 200 };
		
		Map<Integer, SimulatorElement> elems = UCERF3ComparisonAnalysis.loadElements(sol.getRupSet());
		
		List<List<SimulatorEvent>> eventLists = loadUCERF3EtasCatalogs(catalogs, sol, region, elems);
		List<SimulatorEvent> stitched = UCERF3ComparisonAnalysis.stitch(eventLists);
		
		for (int windowLen : windowLens) {
			File outputFile = new File(outputDir, "ucerf3_etas_"+windowLen+"yr.bin");
			SimulatorMomRateVarCalc.writeMomRateTimeSeries(windowLen, stitched, outputFile);
		}
	}

}
