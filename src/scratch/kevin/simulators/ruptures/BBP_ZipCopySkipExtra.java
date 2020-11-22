package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.Enumeration;
import java.util.HashSet;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipArchiveOutputStream;
import org.apache.commons.compress.archivers.zip.ZipFile;

import com.google.common.base.Preconditions;

public class BBP_ZipCopySkipExtra {

	public static void main(String[] args) throws IOException {
		File bbpDir = new File("/data/kevin/bbp/parallel");
		
		File inputDir = new File(bbpDir,
				"2020_05_05-rundir4983_stitched-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-cs500Sites");
		File outputDir = new File(bbpDir,
				"2020_09_03-rundir4983_stitched-all-m6.5-skipYears65000-noHF-vmLA_BASIN_500-cs500Sites");
		int minEventID = 1109191;
		
		File inputFile = new File(inputDir, "results_rotD.zip");
		File outputFile = new File(outputDir, inputFile.getName());
		
		ZipFile inputZip = new ZipFile(inputFile);
		
		ZipArchiveOutputStream out = new ZipArchiveOutputStream(outputFile);
		
		HashSet<Integer> eventIDs = new HashSet<>();
		int entryCount = 0;
		
		Enumeration<? extends ZipArchiveEntry> entries = inputZip.getEntries();
		while (entries.hasMoreElements()) {
			ZipArchiveEntry e = entries.nextElement();
			String name = e.getName();
			Preconditions.checkState(name.startsWith("event_"));
			String subName = name.substring(name.indexOf("_")+1);
			if (name.contains("/"))
				subName = subName.substring(0, subName.indexOf("/"));
			int eventID = Integer.parseInt(subName);
			if (eventID < minEventID)
				continue;
			eventIDs.add(eventID);
			ZipArchiveEntry outEntry = new ZipArchiveEntry(e.getName());
			outEntry.setCompressedSize(e.getCompressedSize());
			outEntry.setCrc(e.getCrc());
			outEntry.setExternalAttributes(e.getExternalAttributes());
			outEntry.setExtra(e.getExtra());
			outEntry.setExtraFields(e.getExtraFields());
			outEntry.setGeneralPurposeBit(e.getGeneralPurposeBit());
			outEntry.setInternalAttributes(e.getInternalAttributes());
			outEntry.setMethod(e.getMethod());
			outEntry.setRawFlag(e.getRawFlag());
			outEntry.setSize(e.getSize());
			out.addRawArchiveEntry(outEntry, inputZip.getRawInputStream(e));
			entryCount++;
		}
		
		System.out.println("Copied "+entryCount+" entries for "+eventIDs.size()+" events");
		
		inputZip.close();
		out.close();
	}

}
