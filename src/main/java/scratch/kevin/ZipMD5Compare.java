package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.io.archive.ArchiveInput;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;

public class ZipMD5Compare {

	public static void main(String[] args) throws IOException {
		File zip1 = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/results_simplified.zip");
		File zip2 = new File("/tmp/results_simplified.zip");
		
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
		System.out.println("Calulating MD5s for Zip1: "+zip1.getAbsolutePath());
		Map<String, String> md5s1 = loadCalcMD5s(zip1, exec);
		System.out.println("Calulating MD5s for Zip2: "+zip2.getAbsolutePath());
		Map<String, String> md5s2 = loadCalcMD5s(zip2, exec);
		
		for (String entry : md5s1.keySet()) {
			String md51 = md5s1.get(entry);
			String md52 = md5s2.get(entry);
			if (md52 == null)
				System.err.println("Zip2 is missing:\t"+entry);
			else if (!md51.equals(md52))
				System.err.println("Zip2 differs for:\t"+entry+"\t("+md51+" != "+md52+")");
		}
		for (String entry : md5s2.keySet())
			if (!md5s1.containsKey(entry))
				System.err.println("Zip1 is missing:\t"+entry);
		
		exec.shutdown();
	}
	
	public static Map<String, String> loadCalcMD5s(File zipFile, ExecutorService exec) throws IOException {
		ArchiveInput input = ArchiveInput.getDefaultInput(zipFile);
		Map<String, String> ret = loadCalcMD5s(input, exec);
		input.close();
		return ret;
	}
	
	public static Map<String, String> loadCalcMD5s(ArchiveInput input, ExecutorService exec) throws IOException {
		ArrayDeque<MessageDigest> mds = new ArrayDeque<>();
		
		List<String> entries = input.entryStream().collect(Collectors.toList());
		List<Future<byte[]>> futures = new ArrayList<>(entries.size());
		
		for (String entry : entries) {
			futures.add(exec.submit(new Callable<byte[]>() {

				@Override
				public byte[] call() throws Exception {
					MessageDigest md;
					synchronized (mds) {
						md = mds.poll();
					}
					if (md == null) {
						try {
							md = MessageDigest.getInstance("MD5");
						} catch (NoSuchAlgorithmException e) {
							throw ExceptionUtils.asRuntimeException(e);
						}
					}

					byte[] buffer = new byte[8192]; // Read in chunks of 8KB
					int bytesRead;
					
					InputStream is;
					synchronized (input) {
						is = input.getInputStream(entry);
					}

					// Read the input stream in chunks and update the digest
					while ((bytesRead = is.read(buffer)) != -1) {
						md.update(buffer, 0, bytesRead);
					}
					
					byte[] ret = md.digest();
					md.reset();
					synchronized (mds) {
						mds.push(md);
					}
					return ret;
				}
			}));
		}
		
		Map<String, String> ret = new HashMap<>(entries.size());
		try {
			for (int i=0; i<entries.size(); i++)
				ret.put(entries.get(i), hashBytesToString(futures.get(i).get()));
		} catch (InterruptedException | ExecutionException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		return ret;
	}
	
	private static String hashBytesToString(byte[] hash) {
		BigInteger bigInt = new BigInteger(1,hash);
		String hashtext = bigInt.toString(16);
		// Now we need to zero pad it if you actually want the full 32 chars.
		while(hashtext.length() < 32 ){
			hashtext = "0"+hashtext;
		}
		return hashtext;
	}

}
