package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.XMLUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.LastEventData;

public class DateLastEventAdder {

	public static void main(String[] args) throws IOException, DocumentException {
		Map<Integer, List<LastEventData>> data = LastEventData.load();
		
		File file;
		if (args.length == 1)
			file = new File(args[0]);
		else
//			file = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/UCERF3_ERF/"
//				+ "cached_FM3_1_dep100.0_depMean_rakeMean.zip");
			file = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/UCERF3_ERF");
		
		if (file.isDirectory()) {
			for (File f : file.listFiles()) {
				if (f.isFile() && f.getName().endsWith(".zip"))
					addIt(f, data);
			}
		} else {
			addIt(file, data);
		}
	}
	
	private static void addIt(File fssFile, Map<Integer, List<LastEventData>> data) throws IOException, DocumentException {
		File tempDir = Files.createTempDir();
		
		System.out.println("Unzipping "+fssFile.getName()+" to "+tempDir.getAbsolutePath());
		FileUtils.unzipFile(fssFile, tempDir);
		
		File fsdFile = new File(tempDir, "fault_sections.xml");
		if (!fsdFile.exists()) {
			System.out.println("Not a sol file, skipping");
			FileUtils.deleteRecursive(tempDir);
			return;
		}
		Preconditions.checkState(fsdFile.exists());
		
		Document doc = XMLUtils.loadDocument(fsdFile);
		Element fsEl = doc.getRootElement().element(FaultSectionPrefData.XML_METADATA_NAME+"List");
		ArrayList<FaultSectionPrefData> faultSectionData = FaultSystemIO.fsDataFromXML(fsEl);
		int withBefore = countLastEventData(faultSectionData);
		LastEventData.populateSubSects(faultSectionData, data);
		int withAfter = countLastEventData(faultSectionData);
		System.out.println("Sects with last event data:");
		System.out.println("\tBefore: "+withBefore+"/"+faultSectionData.size());
		System.out.println("\tAfter: "+withAfter+"/"+faultSectionData.size());
		if (withAfter == withBefore) {
			System.out.println("Identical, skipping");
			FileUtils.deleteRecursive(tempDir);
			return;
		}
		FaultModels fm = null;
		DeformationModels dm = null;
		if (fsEl.attributeValue("faultModName") != null)
			fm = FaultModels.valueOf(fsEl.attributeValue("faultModName"));
		if (fsEl.attributeValue("defModName") != null)
			dm = DeformationModels.valueOf(fsEl.attributeValue("defModName"));
		
		doc = XMLUtils.createDocumentWithRoot();
		Element root = doc.getRootElement();
		FaultSystemIO.fsDataToXML(root, FaultSectionPrefData.XML_METADATA_NAME+"List", fm, dm, faultSectionData);
		XMLUtils.writeDocumentToFile(fsdFile, doc);
		
		List<String> fileNames = Lists.newArrayList();
		for (File file : tempDir.listFiles())
			if (file.isFile())
				fileNames.add(file.getName());
		
		FileUtils.createZipFile(fssFile.getAbsolutePath(), tempDir.getAbsolutePath(), fileNames);
		
		FileUtils.deleteRecursive(tempDir);
	}
	
	private static int countLastEventData(List<FaultSectionPrefData> fsd) {
		int cnt = 0;
		
		for (FaultSectionPrefData sect : fsd)
			if (sect.getDateOfLastEvent() > Long.MIN_VALUE)
				cnt++;
		
		return cnt;
	}

}
