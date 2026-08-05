package scratch.kevin.nshm27.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import org.opensha.sha.util.TectonicRegionType;

import net.mahdilamb.colormap.Colors;

public class NSHM27_PaperPaths {
	
	public static final File PAPER_DIR =new File("/home/kevin/Documents/papers/2027_GNMI_AmSam_ERF/");
	public static final File FIGURES_DIR =new File(PAPER_DIR, "Figures");
	
	public static final File INV_DIR = new File("/home/kevin/OpenSHA/fss_inversions/");
	
	public static final String MODEL_DATE = "2026_07_13";
	
	public static final File AMSAM_SOL_DIR = new File(INV_DIR, MODEL_DATE+"-nshm27-AMSAM-5000samples-lhs_pairwise");
	public static final File GNMI_SOL_DIR = new File(INV_DIR, MODEL_DATE+"-nshm27-GNMI-5000samples-lhs_pairwise");
	
	public static File getSolDir(NSHM27_SeismicityRegions seisReg) {
		return switch (seisReg) {
		case AMSAM:
			yield AMSAM_SOL_DIR;
		case GNMI:
			yield GNMI_SOL_DIR;
		default:
			throw new IllegalArgumentException("Unexpected value: " + seisReg);
		};
	}
	
	public static File getSolFile(NSHM27_SeismicityRegions seisReg) {
		return new File(getSolDir(seisReg), "results_branch_averaged.zip");
	}
	
	public static FaultSystemSolution getSolution(NSHM27_SeismicityRegions seisReg) throws IOException {
		return FaultSystemSolution.load(getSolFile(seisReg));
	}
	
	public static File getInterfaceSolFile(NSHM27_SeismicityRegions seisReg) {
		return new File(getSolDir(seisReg), "results_"+seisReg.name()+"_V1_SUBDUCTION_INTERFACE_branch_averaged.zip");
	}
	
	public static FaultSystemSolution getInterfaceSolution(NSHM27_SeismicityRegions seisReg) throws IOException {
		return FaultSystemSolution.load(getInterfaceSolFile(seisReg));
	}
	
	public static LogicTree<LogicTreeNode> getLogicTree(NSHM27_SeismicityRegions seisReg) throws IOException {
		File ltFile = new File(getSolDir(seisReg), "logic_tree.json");
		return LogicTree.read(ltFile);
	}
	
	public static LogicTree<LogicTreeNode> getAnalysisLogicTree(NSHM27_SeismicityRegions seisReg) throws IOException {
		File ltFile = new File(getSolDir(seisReg), "logic_tree_analysis.json");
		return LogicTree.read(ltFile);
	}
	
	public static TectonicRegionType[] TRTs = {
			TectonicRegionType.SUBDUCTION_INTERFACE,
			TectonicRegionType.SUBDUCTION_SLAB,
			TectonicRegionType.ACTIVE_SHALLOW,
	};
	
	public static Color getColor(TectonicRegionType trt, IncludeBackgroundOption bgType) {
		return switch (trt) {
		case ACTIVE_SHALLOW:
			if (bgType == IncludeBackgroundOption.ONLY)
				yield Colors.tab_lightblue;
			else if (bgType == IncludeBackgroundOption.INCLUDE)
				yield darker(Colors.tab_blue);
			else
				yield Colors.tab_blue;
		case SUBDUCTION_INTERFACE:
			if (bgType == IncludeBackgroundOption.ONLY)
				yield Colors.tab_lightorange;
			else if (bgType == IncludeBackgroundOption.INCLUDE)
				yield darker(Colors.tab_orange);
			else
				yield Colors.tab_orange;
		case SUBDUCTION_SLAB:
			if (bgType == IncludeBackgroundOption.ONLY)
				yield Colors.tab_lightgreen;
			else if (bgType == IncludeBackgroundOption.INCLUDE)
				yield darker(Colors.tab_green);
			else
				yield Colors.tab_green;
		default:
			throw new IllegalArgumentException("Unexpected value: " + trt);
		};
	}
	
	public static Color OBS_RATE_COLOR = Colors.tab_brown.darker();
	
	public static Color darker(Color color) {
		return color.darker();
	}

	public static DecimalFormat oneDF = new DecimalFormat("0.0");
	public static DecimalFormat twoDF = new DecimalFormat("0.00");
	public static DecimalFormat oDF = new DecimalFormat("0.#");

}
