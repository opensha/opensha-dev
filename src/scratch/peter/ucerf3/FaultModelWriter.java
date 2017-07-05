package scratch.peter.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.dom4j.Document;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.util.XMLUtils;
import org.opensha.nshmp2.util.GridUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.FaultSystemIO;

/*
 * Used to write out fault models for local caching
 */
class FaultModelWriter {

	public static void writeFaultModel() {
		try {
			File dir = new File("tmp");
			FaultModels fm = FaultModels.FM3_1;
			ArrayList<FaultSectionPrefData> datas = fm.fetchFaultSections();
			Document doc = XMLUtils.createDocumentWithRoot();
			FaultSystemIO.fsDataToXML(doc.getRootElement(),
				FaultModels.XML_ELEMENT_NAME, fm, null, datas);
			XMLUtils.writeDocumentToFile(new File(dir, fm.getShortName() +
				".xml"), doc);
			System.exit(0);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
//		writeFaultModel();
		writeFaultModelToKML();
//		double[] dd = { 1,2,3,4,5,6,7,8};
//		List<Double> dList = Lists.newArrayList();
//		for (double d : dd) {
//			dList.add(d);
//		}
//		Iterator<Double> it = dList.iterator();
//		while (it.hasNext()) {
//			Double D = it.next();
//			System.out.println(D + " " +D.equals(3.0));
//			if (D.equals(3.0)) it.remove();
//		}
//		System.out.println(dList);
	}
	
	public static void writeFaultModelToKML() {
			FaultModels fm = FaultModels.FM3_2;
			ArrayList<FaultSectionPrefData> datas = fm.fetchFaultSections();
			Map<String, LocationList> traceMap = Maps.newHashMap();
			for (FaultSectionPrefData fspd : datas) {
				traceMap.put(fspd.getSectionName(), fspd.getFaultTrace());
			}
			GridUtils.tracesToKML(traceMap, "UC3_FM32", 
				Color.RED);
			System.exit(0);
	}


}
