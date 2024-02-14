package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CompletableFuture;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureTreeNavigator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;

public class MarkBerkeleyCSVWriter {

	public static void main(String[] args) throws IOException {
		Location loc = new Location(37.8707, -122.2508);
		
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-NSHM23_v3/");
		FaultSystemSolution sol = FaultSystemSolution.load(new File(dir, "results_NSHM23_v3_branch_averaged_gridded.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		File outputDir = new File("/tmp/berkeley_participating_rups");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		boolean doLogicTree = true;
		SolutionLogicTree slt = null;
		if (doLogicTree)
			slt = SolutionLogicTree.load(new File(dir, "results.zip"));
		
		double minDist = Double.POSITIVE_INFINITY;
		FaultSection closestSect = null;
		double traceConsiderDist = 50d;
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			boolean consider = false;
			for (Location traceLoc : sect.getFaultTrace()) {
				if (LocationUtils.horzDistanceFast(traceLoc, loc) < traceConsiderDist) {
					consider = true;
					break;
				}
			}
			if (consider) {
				for (Location traceLoc : sect.getFaultSurface(1d).getEvenlyDiscritizedUpperEdge()) {
					double dist = LocationUtils.horzDistanceFast(traceLoc, loc);
					if (dist < minDist) {
						minDist = dist;
						closestSect = sect;
					}
				}
			}
		}
		Preconditions.checkNotNull(closestSect);
		System.out.println("Closest section: "+closestSect.getSectionId()+". "
				+closestSect.getSectionName()+" ("+(float)minDist+" km away)");
		
		double little_l_Direction = 320d;
		FaultTrace closestTrace = closestSect.getFaultTrace();
		double firstToLastDirection = LocationUtils.azimuth(closestSect.getFaultTrace().first(), closestTrace.last());
		double lastToFirstDirection = LocationUtils.azimuth(closestSect.getFaultTrace().last(), closestTrace.first());
		double diffFirstToLast = angleDiff(little_l_Direction, firstToLastDirection);
		double diffLastToFirst = angleDiff(little_l_Direction, lastToFirstDirection);
		boolean firstTraceLocLittle_l = diffLastToFirst < diffFirstToLast;

		System.out.println("Trace: "+closestTrace);
		System.out.println("First to Last dir: "+(float)firstToLastDirection);
		System.out.println("Last to First dir: "+(float)lastToFirstDirection);
		System.out.println("firstTraceLocLittle_l ? "+firstTraceLocLittle_l);
		FaultTrace resampledTrace = FaultUtils.resampleTrace(closestTrace, 100);
		int closestResampledIndex = -1;
		double closestResampledDist = Double.POSITIVE_INFINITY;
		for (int i=0; i<resampledTrace.size(); i++) {
			double dist = LocationUtils.horzDistanceFast(resampledTrace.get(i), loc);
			if (dist < closestResampledDist) {
				closestResampledDist = dist;
				closestResampledIndex = i;
			}
		}
		double firstSideFract = (double)closestResampledIndex/(double)(resampledTrace.size()-1);
		double little_l_sectDist = firstTraceLocLittle_l ?
				firstSideFract*closestTrace.getTraceLength() : (1d-firstSideFract)*closestTrace.getTraceLength();
		System.out.println("l sectDist = "+(float)little_l_sectDist+" (tot="+(float)closestTrace.getTraceLength()+")");
		
		// claculate lengths
		LengthResult[] lengthResults = new LengthResult[rupSet.getNumRuptures()];
		ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
		List<Integer> rupsForSect = rupSet.getRupturesForSection(closestSect.getSectionId());
		BitSet debugRups = new BitSet(rupSet.getNumRuptures());
		debugRups.set(11433);
		List<Integer> tmpRups = new ArrayList<>(rupsForSect);
		Collections.shuffle(tmpRups, new Random(tmpRups.size()));
		tmpRups = tmpRups.subList(0, 20);
		for (int tmpRup : tmpRups)
			debugRups.set(tmpRup);
		Location closestSectLittleLSide;
		Location closestSectOtherSide;
		if (firstTraceLocLittle_l) {
			closestSectLittleLSide = closestTrace.first();
			closestSectOtherSide = closestTrace.last();
		} else {
			closestSectLittleLSide = closestTrace.last();
			closestSectOtherSide = closestTrace.first();
		}
		DecimalFormat pDF = new DecimalFormat("0.0%");
		for (int rupIndex : rupsForSect) {
			lengthResults[rupIndex] = new LengthResult(closestSect, little_l_sectDist,
					closestSectLittleLSide, closestSectOtherSide, cRups.get(rupIndex));
			if (debugRups.get(rupIndex)) {
				System.out.println("Rup "+rupIndex);
				System.out.println("\tl1="+(float)lengthResults[rupIndex].l1
						+" ("+pDF.format(lengthResults[rupIndex].l1/lengthResults[rupIndex].L)+")");
				System.out.println("\tl2="+(float)lengthResults[rupIndex].l2
						+" ("+pDF.format(lengthResults[rupIndex].l2/lengthResults[rupIndex].L)+")");
				System.out.println("\tL="+(float)lengthResults[rupIndex].L
						+" (rs: "+(float)rupSet.getLengthForRup(rupIndex)+")");
			}
		}
		
		CSVFile<String> csv = buildCSV(sol, closestSect, lengthResults);
		
		csv.writeToFile(new File(outputDir, "participating_ruptures.csv"));
		GeoJSONFaultReader.writeFaultSections(new File(outputDir, "all_sections.geojson"), rupSet.getFaultSectionDataList());
		
		if (doLogicTree) {
			File logicTreeDir = new File(outputDir, "logic_tree_results");
			Preconditions.checkState(logicTreeDir.exists() || logicTreeDir.mkdir());
			
			LogicTree<?> tree = slt.getLogicTree();
			
			CSVFile<String> treeCSV = new CSVFile<>(true);
			List<String> header = new ArrayList<>();
			header.add("File Name");
			header.add("Weight");
			for (LogicTreeLevel<?> level : tree.getLevels())
				header.add(level.getShortName());
			treeCSV.addLine(header);
			
			CompletableFuture<Void> csvFuture = null;
			
			double totWeight = 0d;
			for (LogicTreeBranch<?> branch : tree)
				totWeight += tree.getBranchWeight(branch);
			for (int i=0; i<tree.size(); i++) {
				LogicTreeBranch<?> branch = tree.getBranch(i);
				System.out.println("Branch "+i+": "+branch);
				FaultSystemSolution branchSol = slt.forBranch(branch);
				List<String> line = new ArrayList<>();
				File csvFile = new File(logicTreeDir, branch.buildFileName()+".csv");
				line.add(csvFile.getName());
				line.add((float)(tree.getBranchWeight(branch)/totWeight)+"");
				for (LogicTreeNode node : branch)
					line.add(node.getShortName());
				treeCSV.addLine(line);
				
				if (csvFuture != null)
					csvFuture.join();
				
				FaultSection myClosestSect = closestSect;
				csvFuture = CompletableFuture.runAsync(new Runnable() {
					
					@Override
					public void run() {
						CSVFile<String> csv = buildCSV(branchSol, myClosestSect, lengthResults);
						try {
							csv.writeToFile(csvFile);
						} catch (IOException e) {
							System.exit(1);
						}
					}
				});
			}
			
			if (csvFuture != null)
				csvFuture.join();
			
			treeCSV.writeToFile(new File(logicTreeDir, "logic_tree.csv"));
		}
	}

	public static CSVFile<String> buildCSV(FaultSystemSolution sol, FaultSection closestSect, LengthResult[] lengthResults) {
		CSVFile<String> csv = new CSVFile<>(false);
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		csv.addLine("Rupture Index", "Magnitude", "Average Slip (m)", "Annual Rate", "Average Rake",
				"Average Dip", "l1 (km)", "l2 (km)", "L (km)", "Num Sections", "Section Index 1", "Section Index N");
		
		AveSlipModule slips = rupSet.requireModule(AveSlipModule.class);
		for (int rupIndex : rupSet.getRupturesForSection(closestSect.getSectionId())) {
			List<Integer> sects = rupSet.getSectionsIndicesForRup(rupIndex);
			List<String> line = new ArrayList<>(6+sects.size());
			
			line.add(rupIndex+"");
			line.add((float)rupSet.getMagForRup(rupIndex)+"");
			line.add((float)slips.getAveSlip(rupIndex)+"");
			line.add(sol.getRateForRup(rupIndex)+"");
			line.add((float)rupSet.getAveRakeForRup(rupIndex)+"");
			line.add((float)rupSet.getSurfaceForRupture(rupIndex, 1d).getAveDip()+"");
			line.add((float)lengthResults[rupIndex].l1+"");
			line.add((float)lengthResults[rupIndex].l2+"");
			line.add((float)lengthResults[rupIndex].L+"");
			line.add(sects.size()+"");
			for (int sect : sects)
				line.add(sect+"");
			csv.addLine(line);
		}
		return csv;
	}
	
	private static double angleDiff(double angle1, double angle2) {
		double angleDiff = Math.abs(angle1 - angle2);
		while (angleDiff > 270)
			angleDiff -= 360;
		return Math.abs(angleDiff);
	}
	
	private static class LengthResult {
		final double l1;
		final double l2;
		final double L; // total length
		
		public LengthResult(FaultSection closestSect, double closestSectL1, Location closestSectL1Side,
				Location closestSectL2Side, ClusterRupture rup) {
			RuptureTreeNavigator nav = rup.getTreeNavigator();
			Collection<FaultSection> descendants = nav.getDescendants(closestSect);
			FaultSection predecessor = nav.getPredecessor(closestSect);
			
			double l1DescDist = Double.POSITIVE_INFINITY;
			double l2DescDist = Double.POSITIVE_INFINITY;
			double l1PredDist = Double.POSITIVE_INFINITY;
			double l2PredDist = Double.POSITIVE_INFINITY;
			
			for (FaultSection sect : descendants) {
				l1DescDist = Math.min(l1DescDist,
						Math.min(LocationUtils.horzDistanceFast(closestSectL1Side, sect.getFaultTrace().first()),
								LocationUtils.horzDistanceFast(closestSectL1Side, sect.getFaultTrace().last())));
				l2DescDist = Math.min(l2DescDist,
						Math.min(LocationUtils.horzDistanceFast(closestSectL2Side, sect.getFaultTrace().first()),
								LocationUtils.horzDistanceFast(closestSectL2Side, sect.getFaultTrace().last())));
			}
			if (predecessor != null) {
				l1PredDist = Math.min(l1PredDist,
						Math.min(LocationUtils.horzDistanceFast(closestSectL1Side, predecessor.getFaultTrace().first()),
								LocationUtils.horzDistanceFast(closestSectL1Side, predecessor.getFaultTrace().last())));
				l2PredDist = Math.min(l2PredDist,
						Math.min(LocationUtils.horzDistanceFast(closestSectL2Side, predecessor.getFaultTrace().first()),
								LocationUtils.horzDistanceFast(closestSectL2Side, predecessor.getFaultTrace().last())));
			}
			
			boolean l1DescCloser; // true if little l1 is on the descendant side
			if (Double.isFinite(l1DescDist) && Double.isFinite(l1PredDist)) {
				// section is in the middle of a rupture
				// desc is closer if descendant is closer to l1 than the pred
				l1DescCloser = l1DescDist < l1PredDist;
				if (l1DescCloser)
					Preconditions.checkState(l2PredDist < l2DescDist);
			} else if (Double.isFinite(l1DescDist)) {
				// we only have a descendant side
				// desc is closer if we're closer to l1 than the other end
				l1DescCloser = l1DescDist < l2DescDist;
			} else {
				Preconditions.checkState(Double.isFinite(l1PredDist));
				// we only have a predecessor side
				// desc is closer if we're closer to the other end than l1
				l1DescCloser = l2PredDist < l1PredDist;
			}

			if (l1DescCloser) {
				// l1 is on the descendant (forwards) side
				l1 = calcLength(nav, closestSect, true) + closestSectL1;
				l2 = calcLength(nav, closestSect, false) + (closestSect.getTraceLength()-closestSectL1);
			} else {
				// l1 is on the descendant (forwards) side
				l1 = calcLength(nav, closestSect, false) + (closestSect.getTraceLength()-closestSectL1);
				l2 = calcLength(nav, closestSect, false) + closestSectL1;
			}
			L = l1 + l2;
		}
		
		private double calcLength(RuptureTreeNavigator nav, FaultSection curSect, boolean forwards) {
			List<FaultSection> dests;
			if (forwards)
				dests = List.copyOf(nav.getDescendants(curSect));
			else
				dests = nav.getPredecessor(curSect) == null ? null : List.of(nav.getPredecessor(curSect));
			if (dests == null)
				return 0d;
			
			double sum = 0d;
			for (FaultSection dest : dests) {
				sum += dest.getTraceLength();
				sum += calcLength(nav, dest, forwards);
			}
			return sum;
		}
	}

}
