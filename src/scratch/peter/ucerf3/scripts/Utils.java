package scratch.peter.ucerf3.scripts;

import static com.google.common.base.Charsets.US_ASCII;
import static scratch.UCERF3.enumTreeBranches.DeformationModels.*;
import static scratch.UCERF3.enumTreeBranches.FaultModels.*;
import static scratch.UCERF3.enumTreeBranches.InversionModels.*;
import static scratch.UCERF3.enumTreeBranches.MaxMagOffFault.*;
import static scratch.UCERF3.enumTreeBranches.MomentRateFixes.*;
import static scratch.UCERF3.enumTreeBranches.ScalingRelationships.*;
import static scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels.*;
import static scratch.UCERF3.enumTreeBranches.SpatialSeisPDF.*;
import static scratch.UCERF3.enumTreeBranches.TotalMag5Rate.*;
import static org.opensha.nshmp2.util.Period.*;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.Period;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.enumTreeBranches.MomentRateFixes;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.peter.ucerf3.calc.UC3_CalcUtils;

import com.google.common.base.Charsets;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.io.Files;

/**
 * Add comments here
 * 
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class Utils {

	public static void main(String[] args) throws IOException {
		// generateFullTreeBranchList("UC32tree1440");
		// generateBranchList("UC32-FM-DM-MS-U2-1sec");
		// buildEqnSetWtList();
		curveCheck();
//		mapCheck();
//		invRunCheck(9);
	}

	private static final String LF = IOUtils.LINE_SEPARATOR;

	private static void generateFullTreeBranchList(String fileName) {

		Set<FaultModels> fltModels = EnumSet.of(FM3_1, FM3_2);
		Set<DeformationModels> defModels = EnumSet.of(ABM, GEOLOGIC, NEOKINEMA,
			ZENGBB);
		Set<ScalingRelationships> scalingRel = EnumSet.of(ELLSWORTH_B,
			ELLB_SQRT_LENGTH, HANKS_BAKUN_08, SHAW_CONST_STRESS_DROP,
			SHAW_2009_MOD);
		Set<SlipAlongRuptureModels> slipRup = EnumSet.of(UNIFORM, TAPERED);
		Set<InversionModels> invModels = EnumSet.of(CHAR_CONSTRAINED);
		Set<TotalMag5Rate> totM5rate = EnumSet.of(RATE_6p5, RATE_7p9, RATE_9p6);
		Set<MaxMagOffFault> mMaxOff = EnumSet.of(MAG_7p3, MAG_7p6, MAG_7p9);
		Set<MomentRateFixes> momentFix = EnumSet.of(NONE);
		Set<SpatialSeisPDF> spatialSeis = EnumSet.of(UCERF2, UCERF3);

		buildList(fileName, fltModels, defModels, scalingRel, slipRup,
			invModels, totM5rate, mMaxOff, momentFix, spatialSeis);
	}

	private static void generateBranchList(String fileName) {

		Set<FaultModels> fltModels = EnumSet.of(FM3_1, FM3_2); // FM3_1, FM3_2);
		Set<DeformationModels> defModels = EnumSet.of(ABM, GEOLOGIC, NEOKINEMA,
			ZENGBB); // ABM, GEOLOGIC, NEOKINEMA, ZENGBB);
		Set<ScalingRelationships> scalingRel = EnumSet.of(ELLSWORTH_B,
			ELLB_SQRT_LENGTH, HANKS_BAKUN_08, SHAW_CONST_STRESS_DROP,
			SHAW_2009_MOD); // ELLSWORTH_B, ELLB_SQRT_LENGTH, HANKS_BAKUN_08,
							// SHAW_CONST_STRESS_DROP, SHAW_2009_MOD);
		Set<SlipAlongRuptureModels> slipRup = EnumSet.of(TAPERED); // UNIFORM,
																	// TAPERED);
		Set<InversionModels> invModels = EnumSet.of(CHAR_CONSTRAINED);
		Set<TotalMag5Rate> totM5rate = EnumSet.of(RATE_7p9); // RATE_7p6,
																// RATE_8p7,
																// RATE_10p0);
		Set<MaxMagOffFault> mMaxOff = EnumSet.of(MAG_7p6); // MAG_7p2, MAG_7p6,
															// MAG_8p0);
		Set<MomentRateFixes> momentFix = EnumSet.of(NONE);
		Set<SpatialSeisPDF> spatialSeis = EnumSet.of(UCERF2); // UCERF2,
																// UCERF3);

		buildList(fileName, fltModels, defModels, scalingRel, slipRup,
			invModels, totM5rate, mMaxOff, momentFix, spatialSeis);
	}

	private static void buildList(String fileName, Set<FaultModels> fltModels,
			Set<DeformationModels> defModels,
			Set<ScalingRelationships> scalingRel,
			Set<SlipAlongRuptureModels> slipRup,
			Set<InversionModels> invModels, Set<TotalMag5Rate> totM5rate,
			Set<MaxMagOffFault> mMaxOff, Set<MomentRateFixes> momentFix,
			Set<SpatialSeisPDF> spatialSeis) {

		List<Set<? extends LogicTreeBranchNode<?>>> branchSets = Lists
			.newArrayList();
		branchSets.add(fltModels);
		branchSets.add(defModels);
		branchSets.add(scalingRel);
		branchSets.add(slipRup);
		branchSets.add(invModels);
		branchSets.add(totM5rate);
		branchSets.add(mMaxOff);
		branchSets.add(momentFix);
		branchSets.add(spatialSeis);

		int count = 0;
		Set<List<LogicTreeBranchNode<?>>> branches = Sets
			.cartesianProduct(branchSets);
		try {
			File out = new File("tmp/invSolSets", fileName + ".txt");
			Files.write("", out, US_ASCII);
			for (List<LogicTreeBranchNode<?>> branch : branches) {
				LogicTreeBranch ltb = LogicTreeBranch.fromValues(branch);
				Files.append(ltb.buildFileName() + LF, out, US_ASCII);
				System.out.println((count++) + " " + ltb.buildFileName());
			}
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

	private static void buildEqnSetWtList() throws IOException {
		String path = "tmp/UC33/src/vars/2013_05_09-ucerf3p3-branch-wt-test_COMPOUND_SOL.zip";
		CompoundFaultSystemSolution cfss = UC3_CalcUtils
			.getCompoundSolution(path);
		File out = new File("tmp/invSolSets", "eqnSetWts.txt");
		int count = 0;
		for (LogicTreeBranch branch : cfss.getBranches()) {
			Files.append(branch.buildFileName() + LF, out, US_ASCII);
			System.out.println((count++) + " " + branch.buildFileName());
		}
	}

	// checks to see if curve set is complete and report missing data
	static Set<Period> periods = Sets
		.newHashSet(GM0P00, GM0P20); //, GM1P00, GM4P00);

	private static void curveCheck() throws IOException {
		String srcPath = "tmp/UC33/curves/src/UC33tree13";
		File srcDir = new File(srcPath);
		File[] branchDirs = srcDir.listFiles();
		for (File branchDir : branchDirs) {
			if (!branchDir.isDirectory()) continue;
			for (Period period : periods) {
				File perDir = new File(branchDir, period.name());
				if (!perDir.exists()) {
					System.out.println("Missing dir: " + perDir);
					continue;
				}
				File curveFile = new File(perDir, "curves.csv");
				List<String> lines = Files.readLines(curveFile,
					Charsets.US_ASCII);
				if (lines.size() != 61) {
					System.out.println("Missing curves: " + perDir);
					System.out.println("Curve count: " + lines.size());
				}
			}
		}
	}

	// checks to see if curve set is complete and report missing data
	private static void mapCheck() throws IOException {
		String srcPath = "tmp/UC33/maps/src/UC33";
		File srcDir = new File(srcPath);
		File[] branchDirs = srcDir.listFiles();
		for (File branchDir : branchDirs) {
			if (!branchDir.isDirectory()) continue;

			File gridDir = new File(branchDir, TestGrid.CA_RELM.name());
			if (!gridDir.exists()) {
				System.out.println("Missing dir: " + gridDir);
				continue;
			}
			
			File perDir = new File(gridDir, GM0P00.name());
			if (!perDir.exists()) {
				System.out.println("Missing dir: " + perDir);
				continue;
			}
			
			File curveFile = new File(perDir, "curves.csv");
			List<String> lines = Files.readLines(curveFile, Charsets.US_ASCII);
			if (lines.size() != 7637) {
				System.out.println("Missing curves: " + perDir);
				System.out.println("Curve count: " + lines.size());
			}
			if (!lines.get(0).startsWith("lat")) {
				System.out.println("Malformed curves: " + perDir);
				System.out.println("Curve count: " + lines.size());
			}
		}
	}
	
	
	private static void invRunCheck(int runID) throws IOException {
		Set<Period> periods = EnumSet.of(GM0P00, GM3P00);
		
		String root = "tmp/UC33/curves/src/invRuns/";
		File srcDir = new File(root + "run" + runID);
		File[] branchDirs = srcDir.listFiles();
		
		
		String brPath = "tmp/UC33/curvejobs/branches/tree1440.txt";
		File brFile = new File(brPath);
		List<String> brList = Files.readLines(brFile, US_ASCII);
		Set<String> brAll = Sets.newHashSet(brList);
		Set<String> brRun = Sets.newHashSet();
		Set<String> brErr = Sets.newHashSet();
		
		// loop existing dirs
		// 		check that each dir has GM0P00 and GM3P00
		//		check that each has 41 lines
		
		for (File branchDir : branchDirs) {
			if (!branchDir.isDirectory()) continue;
			String brID = branchDir.getName();
			brRun.add(brID);
			for (Period period : periods) {
				File perDir = new File(branchDir, period.name());
				
				// keep branch if period dir is missing
				if (!perDir.exists()) {
					System.out.println("Missing period dir: " + period);
					System.out.println("Branch: " + brID);
					brErr.add(brID);
					continue;
				}
				
				// keep branch if either curve file is incomplete
				File curveFile = new File(perDir, "curves.csv");
				List<String> lines = Files.readLines(curveFile,
					Charsets.US_ASCII);
				if (lines.size() != 41) {
					System.out.println("Missing curves: " + period + " " + lines.size());
					System.out.println("Branch: " + brID);
					brErr.add(brID);
				}
			}
		}
		
		// create list of branches to rerun: brAll - (brRun - brErr)
		Set<String> rerun = Sets.difference(brAll, Sets.difference(brRun, brErr));
		System.out.println("Total processed: " + brRun.size());
		System.out.println("Errors: " + brErr.size());
		System.out.println("To rerun: " + rerun.size());
		
		File runFix = new File(root, "run" + runID + "fix.txt");
		Files.write("", runFix, US_ASCII);
		for (String br : rerun) {
			Files.append(br + LF, runFix, US_ASCII);
		}
		
		for (String br : brErr) {
			File brDir = new File(srcDir, br);
			FileUtils.deleteDirectory(brDir);
		}
		//		for (String err : brErr) {
//		
//	}

	}


}
