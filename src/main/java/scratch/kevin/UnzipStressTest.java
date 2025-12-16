package scratch.kevin;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.UnknownHostException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;

import org.opensha.commons.util.io.archive.ArchiveInput;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.commons.util.modules.OpenSHA_Module;

import com.google.common.io.Files;

public class UnzipStressTest {

	public static void main(String[] args) {
		if (args.length != 3) {
			System.err.println("USAGE: <file-list.txt> <output-log-file.txt> <trial-count>");
			System.exit(2);
		}
		
		String hostname = null;
		try {
			hostname = java.net.InetAddress.getLocalHost().getHostName();
		} catch (UnknownHostException e) {
			e.printStackTrace();
			System.err.println("Couldn't detect hostname");
			System.exit(2);
		}
		
		List<File> inFiles = null;
		try {
			List<String> lines = Files.readLines(new File(args[0]), Charset.defaultCharset());
			inFiles = new ArrayList<>();
			for (String line : lines) {
				if (line.isBlank())
					continue;
				File inFile = new File(line);
				if (!inFile.isDirectory())
					inFiles.add(inFile);
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("Couldn't read file list");
			System.exit(2);
		}
		
		File logFile = new File(args[1]);
		FileWriter fw = null;
		try {
			fw = new FileWriter(logFile);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		int trials = Integer.parseInt(args[2]);
		
		System.out.println("Will do "+trials+" reads, randomly sampled across "+inFiles.size()+" input files");
		
		ModuleArchive.VERBOSE_DEFAULT = false;
		
		SimpleDateFormat sdf = new SimpleDateFormat();
		
		Random r = new Random(combinedSeed(hostname));
		int numSuccess = 0;
		int numFail = 0;
		for (int i=0; i<trials; i++) {
			int index = r.nextInt(inFiles.size());
			File inFile = inFiles.get(index);
			System.out.println("Trial "+i+"/"+trials+":\t"+inFile.getAbsolutePath());
			try {
				ArchiveInput input = ArchiveInput.getDefaultInput(inFile);
				
				ModuleArchive<OpenSHA_Module> archive = new ModuleArchive<>(input);
				
				List<Callable<OpenSHA_Module>> available = archive.getAvailableModules();
				for (Callable<OpenSHA_Module> call : available)
					call.call();
				
				input.close();
				
				String str = "SUCCESS: "+hostname+" trial "+i+"/"+trials+" reading "+inFile.getAbsolutePath();
				System.out.println("\t"+str);
				fw.write(str);
				fw.write('\n');
				numSuccess++;
			} catch (IOException e) {
				System.err.println("Exception on "+hostname+" trial "+i+"/"+trials+" reading "+inFile.getAbsolutePath()+": "+e.getMessage());
				numFail++;
				try {
					String str = "FAIL: "+hostname+" trial "+i+"/"+trials+" reading "+inFile.getAbsolutePath()+" at "+sdf.format(new Date())+": "+e.getMessage();
					System.out.println("\t"+str);
					fw.write(str);
					fw.write('\n');
				} catch (IOException e1) {
					e1.printStackTrace();
					System.exit(1);
				}
				continue;
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		System.out.println(numSuccess+"/"+trials+" successes");
		System.out.println(numFail+"/"+trials+" failures");
		
		try {
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public static long hostnameSeed(String hostname) {
		try {
			MessageDigest digest = MessageDigest.getInstance("SHA-256");
			byte[] hash = digest.digest(hostname.getBytes(StandardCharsets.UTF_8));
			// fold first 8 bytes into a long
			long value = 0;
			for (int i = 0; i < 8; i++)
				value = (value << 8) | (hash[i] & 0xff);
			return value;
		} catch (NoSuchAlgorithmException e) {
			// fallback
			return hostname.hashCode();
		}
	}

	// splitmix64 mixer for good distribution
	public static long mix(long x) {
		x ^= (x >>> 30) * 0xbf58476d1ce4e5b9L;
		x = (x ^ (x >>> 27)) * 0x94d049bb133111ebL;
		return x ^ (x >>> 31);
	}

	public static long combinedSeed(String hostname) {
		long hostSeed = hostnameSeed(hostname);
		long timeSeed = System.nanoTime();
		return mix(hostSeed ^ timeSeed);
	}

}
