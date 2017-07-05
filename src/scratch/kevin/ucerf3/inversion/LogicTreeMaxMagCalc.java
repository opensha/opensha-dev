package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Joiner;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.SectionClusterList;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class LogicTreeMaxMagCalc {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
//		FaultModels[] faultModels = { FaultModels.FM3_1, FaultModels.FM3_2 };
//		DeformationModels[] defModels = { DeformationModels.ABM, DeformationModels.GEOBOUND, DeformationModels.GEOLOGIC,
//				DeformationModels.GEOLOGIC_PLUS_ABM, DeformationModels.NEOKINEMA, DeformationModels.ZENG };
//		MagAreaRelationships[] magAreas = MagAreaRelationships.values();
//		SlipAlongRuptureModels[] slipAlongs =  {SlipAlongRuptureModels.TAPERED}; // not dependant on this
//		AveSlipForRupModels[] aveSlipModels = { AveSlipForRupModels.ELLSWORTH_B,
//				AveSlipForRupModels.SHAW12_SQRT_LENGTH, AveSlipForRupModels.SHAW_12_CONST_STRESS_DROP,
//				AveSlipForRupModels.SHAW_2009_MOD };
//		
//		InversionModels im = InversionModels.CHAR_CONSTRAINED;
//		
//		double defaultAseis = 0.1;
//		
//		CSVFile<String> csv = new CSVFile<String>(true);
//		
//		csv.addLine("Fault Model", "Deformation Model", "M(A) Relationships", "Dr", "Dsr", "Min Mag", "Mag Mag");
//		
//		for (FaultModels fm : faultModels) {
//			DeformationModels filterBasis = fm.getFilterBasis();
//			DeformationModelFetcher filterBasisFetcher = new DeformationModelFetcher(fm, filterBasis,
//					UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, defaultAseis);
//			LaughTestFilter laughTest = LaughTestFilter.getDefault();
//			SectionClusterList clusters = new SectionClusterList(filterBasisFetcher, laughTest);
//			
//			for (DeformationModels dm : defModels) {
//				if (!dm.isApplicableTo(fm))
//					continue;
//				List<FaultSectionPrefData> faultSectionData = new DeformationModelFetcher(fm, dm,
//						UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, defaultAseis).getSubSectionList();
//				for (MagAreaRelationships ma : magAreas) {
//					for (SlipAlongRuptureModels sal : slipAlongs) {
//						for (AveSlipForRupModels as : aveSlipModels) {
//							InversionFaultSystemRupSet rupSet = new InversionFaultSystemRupSet(
//									clusters, dm, faultSectionData, ma.getMagAreaRelationships(), im, sal, as,
//									8.7, 7.6, false, SpatialSeisPDF.UCERF3);
//							double[] mags = rupSet.getMagForAllRups();
//							double max = StatUtils.max(mags);
//							double min = StatUtils.min(mags);
//							csv.addLine(fm.toString(), dm.toString(), ma.toString(), as.toString(), sal.toString(), min+"", max+"");
//							System.out.println(Joiner.on(", ").join(csv.getLine(csv.getNumRows()-1)));
//						}
//					}
//				}
//			}
//		}
//		csv.writeToFile(new File("/tmp/logic_tree_mags.csv"));
//		csv.writeToTabSeparatedFile(new File("/tmp/logic_tree_mags.txt"), 1);
	}

}
