package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.util.DataUtils;
import org.opensha.refFaultParamDb.vo.DeformationModel;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.BranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;

public class BranchAveragedSubSeismoTest {

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		double totWeight = 0d;
		
		FaultModels fm = FaultModels.FM3_1;
		File invDir = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions");
		BranchWeightProvider weightProv = new APrioriBranchWeightProvider();
		boolean only20 = true;
		HashSet<String> done20 = new HashSet<String>();
		
		EvenlyDiscretizedFunc[] subSeismos = null;
		
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File(invDir, "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		
		InversionFaultSystemSolution baSol = FaultSystemIO.loadInvSol(
				new File(invDir, "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		int done = 0;
		for (LogicTreeBranch branch : cfss.getBranches()) {
			if (branch.getValue(FaultModels.class) != fm)
				continue;
			if (only20) {
				String str = branch.getValue(DeformationModels.class).name()
						+"_"+branch.getValue(ScalingRelationships.class).name();
				if (done20.contains(str))
					continue;
				else
					done20.add(str);
			}
			InversionFaultSystemSolution sol = cfss.getSolution(branch);
			double weight = weightProv.getWeight(branch);
			
			if (subSeismos == null)
				subSeismos = new EvenlyDiscretizedFunc[sol.getRupSet().getNumSections()];
			
			List<GutenbergRichterMagFreqDist> subSeismoList = sol.getFinalSubSeismoOnFaultMFD_List();
			for (int i=0; i<subSeismoList.size(); i++) {
				GutenbergRichterMagFreqDist subSeismo = subSeismoList.get(i);
				if (subSeismos[i] == null)
					subSeismos[i] = new EvenlyDiscretizedFunc(
							subSeismo.getMinX(), subSeismo.size(), subSeismo.getDelta());
				else
					Preconditions.checkState(subSeismo.size() == subSeismos[i].size()
							&& (float)subSeismo.getMinX() == (float)subSeismos[i].getMinX());
				
				for (int j=0; j<subSeismo.size(); j++)
					subSeismos[i].set(j, subSeismos[i].getY(j)+weight*subSeismo.getY(j));
			}
			
			done++;
			totWeight += weight;
			
			if (done % 10 == 0)
				System.out.println("Done with "+done);
		}
		
		// now scale by tot weight
		for (int i=0; i<subSeismos.length; i++)
			subSeismos[i].scale(1d/totWeight);
		
		List<GutenbergRichterMagFreqDist> baFuncs = baSol.getFinalSubSeismoOnFaultMFD_List();
		// now check differences
		double totMaxAbsDiff = 0;
		double totMaxPDiff = 0;
		
		int biggestAbsIndex = -1;
		int biggestPercentIndex = -1;
		for (int i=0; i<subSeismos.length; i++) {
			EvenlyDiscretizedFunc meanFunc = subSeismos[i];
			EvenlyDiscretizedFunc baFunc = baFuncs.get(i);
			
			double maxAbsDiff = 0;
			double maxPDiff = 0;
			
			for (int j=0; j<meanFunc.size(); j++) {
				double meanY = meanFunc.getY(j);
				double baY = baFunc.getY(j);
				
				double absDiff = Math.abs(meanY - baY);
				double pDiff = DataUtils.getPercentDiff(baY, meanY);
				
				if (absDiff > maxAbsDiff)
					maxAbsDiff = absDiff;
				if (pDiff > maxPDiff)
					maxPDiff = pDiff;
			}
			
			if (maxAbsDiff > totMaxAbsDiff) {
				totMaxAbsDiff = maxAbsDiff;
				biggestAbsIndex = i;
			}
			if (maxPDiff > totMaxPDiff) {
				totMaxPDiff = maxPDiff;
				biggestPercentIndex = i;
			}
			
			System.out.println(baSol.getRupSet().getFaultSectionData(i).getName()
					+": max_diff: "+maxAbsDiff+"\t"+maxPDiff+" %");
		}
		
		System.out.println("TOTAL: max_diff: "+totMaxAbsDiff+"\t"+totMaxPDiff+" %");
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		funcs.add(subSeismos[biggestAbsIndex]);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		funcs.add(baFuncs.get(biggestAbsIndex));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		new GraphWindow(funcs, baSol.getRupSet().getFaultSectionData(biggestAbsIndex).getName()
				+" MFDs, biggest abs discrep="+totMaxAbsDiff, chars);
		
		funcs = Lists.newArrayList();
		chars = Lists.newArrayList();
		funcs.add(subSeismos[biggestPercentIndex]);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		funcs.add(baFuncs.get(biggestPercentIndex));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		new GraphWindow(funcs, baSol.getRupSet().getFaultSectionData(biggestPercentIndex).getName()
				+" MFDs, biggest percent discrep="+totMaxPDiff, chars);
	}

}
