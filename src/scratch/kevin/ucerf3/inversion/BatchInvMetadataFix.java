package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.Map;
import java.util.concurrent.ConcurrentMap;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.apache.commons.io.IOUtils;
import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.opensha.commons.util.XMLUtils;

import scratch.UCERF3.inversion.BatchPlotGen;
import scratch.UCERF3.inversion.InversionConfiguration;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.VariableLogicTreeBranch;

import com.google.common.base.Preconditions;
import com.google.common.collect.Maps;

public class BatchInvMetadataFix {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws ZipException 
	 * @throws DocumentException 
	 */
	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		File runDir = new File(args[0]);
		File compoundDir = new File(args[1]);
		
		int cnt = 0;
		
		Map<VariableLogicTreeBranch, Map<String, Double>> misfitsMap = Maps
				.newConcurrentMap();
		Map<VariableLogicTreeBranch, Map<String, Double>> energiesMap = Maps
				.newConcurrentMap();
		
		File[] runDirFiles = runDir.listFiles();
		for (File runSubDir : runDirFiles) {
			if (!runSubDir.isDirectory())
				continue;
			String name = runSubDir.getName();
			if (!name.endsWith("_run0"))
				continue;
			VariableLogicTreeBranch branch = VariableLogicTreeBranch.fromFileName(name);
			Preconditions.checkState(branch.isFullySpecified());
			
			File solFile = new File(runSubDir, name+"_sol.zip");
			Preconditions.checkState(solFile.exists());
			ZipFile zip = new ZipFile(solFile);
			ZipEntry invEntry = zip.getEntry("inv_sol_metadata.xml");
			File outputFile = new File(compoundDir, branch.buildFileName()+"_inv_sol_metadata.xml");
			InputStream is = zip.getInputStream(invEntry);
			FileWriter fw = new FileWriter(outputFile);
			IOUtils.copy(is, fw);
			fw.close();
			
			// load misfits
			Document invDoc = XMLUtils.loadDocument(zip.getInputStream(invEntry));
			Element invRoot = invDoc.getRootElement().element("InversionFaultSystemSolution");
			
			Element confEl = invRoot.element(InversionConfiguration.XML_METADATA_NAME);
			InversionConfiguration conf = InversionConfiguration.fromXMLMetadata(confEl);
			
			Element energiesEl = invRoot.element("Energies");
			Map<String, Double> energies = Maps.newHashMap();
			for (Element energyEl : XMLUtils.getSubElementsList(energiesEl)) {
				String type = energyEl.attributeValue("type");
				double value = Double.parseDouble(energyEl.attributeValue("value"));
				energies.put(type, value);
			}
			energiesMap.put(branch, energies);
			Map<String, Double> misfits = InversionFaultSystemSolution.getMisfits(energies, conf);
			misfitsMap.put(branch, misfits);
			
			cnt++;
		}
		System.out.println("Wrote "+cnt+" xml files");
		
		// write energy/misfit CSV
		BatchPlotGen.writeMisfitsCSV(runDir, runDir.getName()+"_energies", energiesMap);
		BatchPlotGen.writeMisfitsCSV(runDir, misfitsMap);
	}

}
