package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.Enumeration;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipArchiveOutputStream;
import org.apache.commons.compress.archivers.zip.ZipFile;

public class ResultsZipRebuildFromAvgGrid {

	public static void main(String[] args) throws IOException {
		File dir = new File("/project/scec_608/kmilner/nshm23/batch_inversions/"
				+ "2022_12_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		File inFile = new File(dir, "results_avg_gridded.zip");
		File outFile = new File(dir, "results.zip");
		
		ZipFile zip = new ZipFile(inFile);
		ZipArchiveOutputStream zout = new ZipArchiveOutputStream(outFile);

		Enumeration<? extends ZipArchiveEntry> entries = zip.getEntries();
		while (entries.hasMoreElements()) {
			ZipArchiveEntry sourceEntry = entries.nextElement();
			String name = sourceEntry.getName();
			if (name.contains("grid_"))
				continue;
			ZipArchiveEntry outEntry = new ZipArchiveEntry(sourceEntry.getName());
			copyEntry(zip, sourceEntry, zout, outEntry);
		}
		
		zip.close();
		zout.close();
	}

	private static void copyEntry(ZipFile sourceZip, ZipArchiveEntry sourceEntry,
			ZipArchiveOutputStream out, ZipArchiveEntry outEntry) throws IOException {
		System.out.println("AsyncLogicTree: copying to zip file: "+outEntry.getName());
		outEntry.setCompressedSize(sourceEntry.getCompressedSize());
		outEntry.setCrc(sourceEntry.getCrc());
		outEntry.setExternalAttributes(sourceEntry.getExternalAttributes());
		outEntry.setExtra(sourceEntry.getExtra());
		outEntry.setExtraFields(sourceEntry.getExtraFields());
		outEntry.setGeneralPurposeBit(sourceEntry.getGeneralPurposeBit());
		outEntry.setInternalAttributes(sourceEntry.getInternalAttributes());
		outEntry.setMethod(sourceEntry.getMethod());
		outEntry.setRawFlag(sourceEntry.getRawFlag());
		outEntry.setSize(sourceEntry.getSize());
		out.addRawArchiveEntry(outEntry, sourceZip.getRawInputStream(sourceEntry));
	}

}
