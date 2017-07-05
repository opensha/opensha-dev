package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.apache.pdfbox.exceptions.COSVisitorException;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.dom4j.DocumentException;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;
import scratch.UCERF3.utils.FaultSystemIO;

public class BulkMagNumPlotGen {

	public static void main(String[] args) throws IOException, DocumentException {
		// assumes each subdirectory contains results and is named by the fss rup index
//		File inputDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_interns/results");
//		File inputDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_interns/results_orig");
		File[] inputDirs = {
				new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_interns/results_orig"),
				new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_interns/results")
		};
		
		File outputDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_interns/plots");
		if (!outputDir.exists())
			outputDir.mkdir();
		
		FaultSystemSolution sol = FaultSystemIO.loadSol(new File(
				"/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		List<File> subDirs = Lists.newArrayList();
		for (File inputDir : inputDirs)
			subDirs.addAll(Lists.newArrayList(inputDir.listFiles()));
		Collections.sort(subDirs, new FileNameComparator());
		
		List<File> pdfFiles =Lists.newArrayList();
		
		for (File dir : subDirs) {
			if (!dir.isDirectory())
				continue;
			String name = dir.getName();
			System.out.println("Processing: "+name);
			int id = Integer.parseInt(name);
			
			double mainshockMag = sol.getRupSet().getMagForRup(id);
			List<String> parentNames = Lists.newArrayList();
			for (FaultSectionPrefData sect : sol.getRupSet().getFaultSectionDataForRupture(id)) {
				String parentName = sect.getParentSectionName();
				if (parentNames.isEmpty() || !parentNames.get(parentNames.size()-1).equals(parentName))
					parentNames.add(parentName);
			}
			String rupName = "M"+(float)mainshockMag+" on "+Joiner.on(", ").join(parentNames);
			
			List<List<ETAS_EqkRupture>> catalogs = Lists.newArrayList();
			for (File subdir : dir.listFiles()) {
				if (!subdir.getName().startsWith("sim_") || !subdir.isDirectory())
					continue;
				if (!MPJ_ETAS_Simulator.isAlreadyDone(subdir))
					continue;
				File catalogFile = new File(subdir, "simulatedEvents.txt");
				catalogs.add(ETAS_CatalogIO.loadCatalog(catalogFile));
			}
			System.out.println("Loaded "+catalogs.size()+" catalogs");
			
			File outputFile = new File(outputDir, name+"_mag_vs_num.png");
			ETAS_SimAnalysisTools.plotMaxMagVsNumAftershocks(catalogs, 0, mainshockMag, outputFile, rupName);
			outputFile = new File(outputDir, name+"_mag_vs_num.pdf");
			ETAS_SimAnalysisTools.plotMaxMagVsNumAftershocks(catalogs, 0, mainshockMag, outputFile, rupName);
			
			pdfFiles.add(outputFile);
		}
		
		PDDocument document = new PDDocument();
		List<PDDocument> subDocs = Lists.newArrayList();
		for (File pdfFile : pdfFiles) {
			PDDocument part = PDDocument.load(pdfFile);
			List<PDPage> list = part.getDocumentCatalog().getAllPages();
			document.addPage(list.get(0));
			subDocs.add(part);
		}
		try {
			document.save(new File(outputDir, "mag_vs_num_combined.pdf"));
		} catch (COSVisitorException e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}
		document.close();
		for (PDDocument doc : subDocs)
			doc.close();
	}

}
