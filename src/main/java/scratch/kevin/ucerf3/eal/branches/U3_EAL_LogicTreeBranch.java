package scratch.kevin.ucerf3.eal.branches;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;

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
		super(build(tiBranch, probModel, gmm, gmmEpi, vs30));
		this.tiBranch = tiBranch;
		this.fssIndexBinFile = fssIndexBinFile;
		this.griddedBinFile = griddedBinFile;
		this.tractDir = tractDir;
	}
	
	private static List<LogicTreeBranchNode<?>> build(U3LogicTreeBranch tiBranch, U3_EAL_ProbModels probModel, U3_EAL_GMMs gmm,
			U3_EAL_GMM_Epistemic gmmEpi, U3_EAL_Vs30Model vs30) {
		List<LogicTreeBranchNode<?>> branches = new ArrayList<>();
		for (LogicTreeBranchNode<?> node : tiBranch)
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
