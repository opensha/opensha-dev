package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.dom4j.Document;
import org.opensha.commons.util.XMLUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.U3FaultSystemIO;

public class FaultModelCacher {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File("/tmp");
//		File dir = new File("D:\\Documents\\temp");
		for (FaultModels fm : FaultModels.values()) {
			ArrayList<? extends FaultSection> datas = fm.fetchFaultSections(true);
			Document doc = XMLUtils.createDocumentWithRoot();
			U3FaultSystemIO.fsDataToXML(doc.getRootElement(), FaultModels.XML_ELEMENT_NAME, fm, null, datas);
			XMLUtils.writeDocumentToFile(new File(dir, fm.getShortName()+".xml"), doc);
		}
		System.exit(0);
	}

}
