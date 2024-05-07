package scratch.kevin.ucerf3.eal.branches;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.logicTree.LogicTreeLevel;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.logicTree.U3LogicTreeBranchNode;

public class U3_EAL_LogicTreeBranch extends U3LogicTreeBranch {
	
	private U3LogicTreeBranch tiBranch;
	private File fssIndexBinFile;
	private File griddedBinFile;
	private File tractDir;

	public U3_EAL_LogicTreeBranch(U3LogicTreeBranch tiBranch, U3_EAL_ProbModels probModel, U3_EAL_GMMs gmm,
			U3_EAL_GMM_Epistemic gmmEpi, U3_EAL_Vs30Model vs30) {
		this(tiBranch, probModel, gmm, gmmEpi, vs30, null, null, null);
	}

	public U3_EAL_LogicTreeBranch(U3LogicTreeBranch tiBranch, U3_EAL_ProbModels probModel, U3_EAL_GMMs gmm,
			U3_EAL_GMM_Epistemic gmmEpi, U3_EAL_Vs30Model vs30, File fssIndexBinFile, File griddedBinFile, File tractDir) {
		super(getLogicTreeLevels(), build(tiBranch, probModel, gmm, gmmEpi, vs30));
		this.tiBranch = tiBranch;
		this.fssIndexBinFile = fssIndexBinFile;
		this.griddedBinFile = griddedBinFile;
		this.tractDir = tractDir;
	}
	
	private static List<LogicTreeLevel<? extends U3LogicTreeBranchNode<?>>> levels;
	
	public static synchronized List<LogicTreeLevel<? extends U3LogicTreeBranchNode<?>>> getLogicTreeLevels() {
		if (levels == null) {
			levels = new ArrayList<>(U3LogicTreeBranch.getLogicTreeLevels());
			List<Class<? extends U3LogicTreeBranchNode<?>>> ealClasses = new ArrayList<>();
			ealClasses.add(U3_EAL_ProbModels.class);
			ealClasses.add(U3_EAL_GMMs.class);
			ealClasses.add(U3_EAL_GMM_Epistemic.class);
			ealClasses.add(U3_EAL_Vs30Model.class);
			for (Class<? extends U3LogicTreeBranchNode<?>> clazz : ealClasses) {
				U3LogicTreeBranchNode<?> value0 = clazz.getEnumConstants()[0];
				LogicTreeLevel<U3LogicTreeBranchNode<?>> level = LogicTreeLevel.forEnumUnchecked(
						value0, value0.getBranchLevelName(), value0.getShortBranchLevelName());
				levels.add(level);
			}
		}
		
		return levels;
	}
	
	private static List<U3LogicTreeBranchNode<?>> build(U3LogicTreeBranch tiBranch, U3_EAL_ProbModels probModel, U3_EAL_GMMs gmm,
			U3_EAL_GMM_Epistemic gmmEpi, U3_EAL_Vs30Model vs30) {
		List<U3LogicTreeBranchNode<?>> branches = new ArrayList<>();
		for (U3LogicTreeBranchNode<?> node : tiBranch)
			branches.add(node);
		branches.add(probModel);
		branches.add(gmm);
		branches.add(gmmEpi);
		branches.add(vs30);
		return branches;
	}
	
	public U3LogicTreeBranch getTIBranch() {
		return tiBranch;
	}

	public File getFSSIndexedBinFile() {
		return fssIndexBinFile;
	}

	public File getGriddedBinFile() {
		return griddedBinFile;
	}

	public File getTractDir() {
		return tractDir;
	}
	
	public static void main(String[] args) {
		U3LogicTreeBranch tiBranch = U3LogicTreeBranch.DEFAULT;
		System.out.println("TI Weight: "+tiBranch.getAprioriBranchWt());
		InversionModels im = tiBranch.getValue(InversionModels.class);
		for (U3_EAL_ProbModels probModel : U3_EAL_ProbModels.values()) {
			System.out.println(probModel.name()+": "+probModel.getRelativeWeight(im));
			for (U3_EAL_GMMs gmm : U3_EAL_GMMs.values()) {
				System.out.println(gmm.name()+": "+gmm.getRelativeWeight(im));
				for (U3_EAL_GMM_Epistemic gmmEpi : U3_EAL_GMM_Epistemic.values()) {
					System.out.println(gmmEpi.name()+": "+gmmEpi.getRelativeWeight(im));
					for (U3_EAL_Vs30Model vs30 : U3_EAL_Vs30Model.values()) {
						System.out.println(vs30.name()+": "+vs30.getRelativeWeight(im));
						U3_EAL_LogicTreeBranch branch = new U3_EAL_LogicTreeBranch(DEFAULT, U3_EAL_ProbModels.POISSON, gmm, gmmEpi, vs30);
						double weight = branch.getAprioriBranchWt();
						System.out.println(branch+"\t\t"+weight);
						Preconditions.checkState(Double.isFinite(weight), "Bad weight: %s", weight);
					}
				}
			}
		}
		
	}

}
