package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.util.modules.ArchiveInput;
import org.opensha.commons.util.modules.ArchiveOutput;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import scratch.kevin.ZipMD5Compare;

public class ParallelZipTests {

	public static void main(String[] args) throws IOException {
		File inFile = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/results.zip");
//		File inFile = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/results_gridded_branches_simplified.zip");
		ArchiveInput input = new ArchiveInput.ZipFileInput(inFile); // don't use apache, we want to de/recompress for this test
//		ModuleArchiveOutput output = new ModuleArchiveOutput.ZipFileOutput(new File("/tmp/benchmark_output-java_zip.zip"));
//		ModuleArchiveOutput output = new ModuleArchiveOutput.ApacheZipFileOutput(new File("/tmp/benchmark_output-apache_zip.zip"));
		ArchiveOutput output = new ArchiveOutput.ParallelZipFileOutput(new File("/tmp/benchmark_output-apache_parallel_zip.zip"), 32, true);
		
		System.out.println("Pre-reading and calculating input MD5s");
		ExecutorService exec = Executors.newFixedThreadPool(8);
		Map<String, String> inputHashes = ZipMD5Compare.loadCalcMD5s(input, exec);
		
		Stopwatch watch = Stopwatch.createStarted();
		if (output instanceof ArchiveOutput.ParallelZipFileOutput)
			((ArchiveOutput.ParallelZipFileOutput)output).setTrackBlockingTimes(true);
		List<String> entries = input.entryStream().toList();
		int numDone = 0;
		for (String entry : entries) {
			System.out.println("Processing "+entry);
			output.transferFrom(input.getInputStream(entry), entry);
			numDone++;
			double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
			double secsEach = secs/(double)numDone;
			double writePerSec = (double)numDone/secs;
			String str = numDone+"/"+entries.size()+" done in "+(float)secs+" s; "+(float)secsEach+" s/write; "+(float)writePerSec+" write/s";
			if (output instanceof ArchiveOutput.ParallelZipFileOutput)
				str += "\t"+((ArchiveOutput.ParallelZipFileOutput)output).getBlockingTimeStats();
			System.out.println(str);
		}
		
		output.close();
		input.close();
		watch.stop();
		ArchiveInput input2 = output.getCompletedInput();
		System.out.println("Calculating output MD5s");
		Map<String, String> outputHashes = ZipMD5Compare.loadCalcMD5s(input2, exec);
		exec.shutdown();
		for (String entry : entries) {
			Preconditions.checkState(outputHashes.containsKey(entry));
			String inputHash = inputHashes.get(entry);
			String outputHash = outputHashes.get(entry);
			Preconditions.checkState(inputHash.equals(outputHash),
					"Bad hash for %s: %s != %s", entry, inputHash, outputHash);
		}
		System.out.println("Validated "+entries.size()+" entries!");
		input2.close();
	}

}
