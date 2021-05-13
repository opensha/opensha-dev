package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.util.Enumeration;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipArchiveOutputStream;
import org.apache.commons.compress.archivers.zip.ZipFile;

public class ZipEntryBatchRename {

	public static void main(String[] args) throws IOException {
		File dir = new File("/data/kevin/bbp/parallel/2019_01_17-rundir2585-rotatedRups-m6p6_vert_ss_surface-"
				+ "50.0km-36srcAz-4siteSrcAz-100rups-skipYears5000-noHF-csLASites");
		File infile = new File(dir, "results_rotD.zip");
		File outfile = new File(dir, "results_rotD_fix.zip");
		
		ZipFile in = new ZipFile(infile);
		ZipArchiveOutputStream out = new ZipArchiveOutputStream(outfile);
		
		Enumeration<? extends ZipArchiveEntry> entries = in.getEntries();
		while (entries.hasMoreElements()) {
			ZipArchiveEntry e = entries.nextElement();
			ZipArchiveEntry outEntry = new ZipArchiveEntry("m6p6_vert_ss_surface_"+e.getName());
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
			out.addRawArchiveEntry(outEntry, in.getRawInputStream(e));
		}
		in.close();
		out.close();
	}

}
