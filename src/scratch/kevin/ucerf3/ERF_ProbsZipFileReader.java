package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import com.google.common.base.Preconditions;

import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.MatrixIO;

public class ERF_ProbsZipFileReader {
	
	private ZipFile zip;
	
	public ERF_ProbsZipFileReader(File zipFile) throws ZipException, IOException {
		this(new ZipFile(zipFile));
	}
	
	public ERF_ProbsZipFileReader(ZipFile zip) {
		this.zip = zip;
	}
	
	/**
	 * Returns probabilities for each rupture, organized by fault system solution rupture index
	 * @param branch
	 * @return
	 * @throws IOException 
	 */
	public double[] getProbabilities(LogicTreeBranch branch) throws IOException {
		String eName = branch.buildFileName()+".bin";
		ZipEntry probsEntry = zip.getEntry(eName);
		Preconditions.checkNotNull(probsEntry, "Entry not found in zip: "+eName);
		
		double[] probs = MatrixIO.doubleArrayFromInputStream(
				zip.getInputStream(probsEntry), probsEntry.getSize());
		
		return probs;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
