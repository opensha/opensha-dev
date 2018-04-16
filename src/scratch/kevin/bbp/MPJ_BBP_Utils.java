package scratch.kevin.bbp;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Enumeration;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipArchiveOutputStream;
import org.apache.commons.compress.archivers.zip.ZipFile;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import edu.usc.kmilner.mpj.taskDispatch.AsyncPostBatchHook;

public class MPJ_BBP_Utils {
	
	public static abstract class MasterZipHook extends AsyncPostBatchHook {
		
		private File zipFile;
		private File rdZipFile;
		
		private ZipArchiveOutputStream out;
		private ZipArchiveOutputStream outRD;
		private byte[] buffer = new byte[1048576];

		public MasterZipHook(File zipFile, File rdZipFile) {
			super(1);
			this.zipFile = zipFile;
			this.rdZipFile = rdZipFile;
		}

		@Override
		protected synchronized void batchProcessedAsync(int[] batch, int processIndex) {
			debug("running async post-batch hook for process "+processIndex+". "+getCountsString());
			debug("async post-batch extimates: "+getRatesString());
			try {
				if (out == null && zipFile != null) {
					if (zipFile.exists())
						Files.move(zipFile, new File(zipFile.getAbsolutePath()+".prev"));
					out = new ZipArchiveOutputStream(new BufferedOutputStream(
							new FileOutputStream(zipFile), buffer.length*4));
				}
				if (outRD == null && rdZipFile != null) {
					if (rdZipFile.exists())
						Files.move(rdZipFile, new File(rdZipFile.getAbsolutePath()+".prev"));
					outRD = new ZipArchiveOutputStream(new BufferedOutputStream(
							new FileOutputStream(rdZipFile), buffer.length*4));
				}
				Preconditions.checkState(out != null || outRD != null);
				for (int index : batch) {
					File subZipFile = getSimZipFile(index);
					String subZipName = subZipFile.getName();
					Preconditions.checkState(subZipFile.exists() && subZipName.endsWith(".zip"));
					
					String simDirName = subZipName.substring(0, subZipName.indexOf(".zip"))+"/";
					ZipArchiveEntry dirEntry = new ZipArchiveEntry(simDirName);
					if (out != null) {
						out.putArchiveEntry(dirEntry);
						out.closeArchiveEntry();
					}
					if (outRD != null) {
						outRD.putArchiveEntry(dirEntry);
						outRD.closeArchiveEntry();
					}
					
					ZipFile sub;
					try {
						sub = new ZipFile(subZipFile);
					} catch (Exception e1) {
						debug("Error with "+subZipFile.getAbsolutePath()+": "+e1.getMessage());
						File subDir = new File(subZipFile.getParentFile(), simDirName);
						if (subDir.exists()) {
							debug("Re-zipping "+simDirName+" from directory");
							FileUtils.createZipFile(subZipFile, subDir, true);
							FileUtils.deleteRecursive(subDir);
							sub = new ZipFile(subZipFile);
						} else {
							throw e1;
						}
					}
					Enumeration<? extends ZipArchiveEntry> entries = sub.getEntries();
					while (entries.hasMoreElements()) {
						ZipArchiveEntry e = entries.nextElement();
						boolean rd = e.getName().endsWith(".rd50") || e.getName().endsWith(".rd100");
						if (out == null && !rd)
							continue;
						ZipArchiveEntry outEntry = new ZipArchiveEntry(simDirName+e.getName());
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
						if (out != null)
							out.addRawArchiveEntry(outEntry, sub.getRawInputStream(e));
						if (rd && outRD != null)
							outRD.addRawArchiveEntry(outEntry, sub.getRawInputStream(e));
					}
					sub.close();
				}
			} catch (Exception e) {
				e.printStackTrace();
				abortAndExit(2);
			}
			debug("done running async post-batch hook for process "+processIndex+". "+getCountsString());
		}

		@Override
		public void shutdown() {
			super.shutdown();
			try {
				if (out != null)
					out.close();
				if (outRD != null)
					outRD.close();
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		
		protected abstract void debug(String message);
		
		protected abstract void abortAndExit(int status);
		
		protected abstract File getSimZipFile(int index);
		
	}
	
	public static Options addCommonOptions(Options ops, boolean includeSitesFile, boolean includeSourceFile,
			boolean includeSRFFile, boolean srfRequired) {
		Option vmOp = new Option("vm", "vm", true, "Velocity model");
		vmOp.setRequired(true);
		ops.addOption(vmOp);
		
		Option methodOp = new Option("m", "method", true, "BBP method");
		methodOp.setRequired(true);
		ops.addOption(methodOp);
		
		if (includeSourceFile) {
			Option srcFile = new Option("src", "src-file", true, "Source file");
			srcFile.setRequired(true);
			ops.addOption(srcFile);
		}
		
		if (includeSRFFile) {
			Option srfFile = new Option("srf", "srf-file", true, "SRF file");
			srfFile.setRequired(srfRequired);
			ops.addOption(srfFile);
		}
		
		if (includeSitesFile) {
			Option sitesFile = new Option("sites", "sites-file", true, "Sites file");
			sitesFile.setRequired(true);
			ops.addOption(sitesFile);
		}
		
		Option outputDir = new Option("o", "output-dir", true, "Output dir");
		outputDir.setRequired(true);
		ops.addOption(outputDir);
		
		Option noHF = new Option("nhf", "no-hf", false, "Flag to disable high-frequency");
		noHF.setRequired(false);
		ops.addOption(noHF);
		
		Option dataDir = new Option("data", "bbp-data-dir", true, "Path to bbp_data dir");
		dataDir.setRequired(false);
		ops.addOption(dataDir);
		
		Option gfDir = new Option("gf", "bbp-gf-dir", true, "Path to bbp_gf dir");
		gfDir.setRequired(false);
		ops.addOption(gfDir);
		
		return ops;
	}
	
	public static void waitOnDir(File dir, int maxRetries, long sleepMillis) {
		int retry = 0;
		while (!(dir.exists() || dir.mkdir())) {
			try {
				Thread.sleep(sleepMillis);
			} catch (InterruptedException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			if (retry++ > maxRetries)
				throw new IllegalStateException("Directory doesn't exist and couldn't be created after "
						+maxRetries+" retries: "+dir.getAbsolutePath());
		}
	}

}
