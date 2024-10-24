package scratch.kevin;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.util.ComparablePairing;

import com.google.common.base.Preconditions;

public class ZipListCompressedDirSizes {

	public static void main(String[] args) throws IOException {
//		ZipFile zip = new ZipFile("/home/kevin/workspace/scec_vdo_vtk/lib/opensha.jar");
		ZipFile zip = new ZipFile("/home/kevin/workspace/opensha/build/libs/opensha-all.jar");
		
		double minSizeToExpand = 1d;
		double minSizeToPrint = 0.1d;
		
		Enumeration<? extends ZipEntry> entries = zip.entries();
		Map<String, Long> pathSizes = new HashMap<>();
		Map<String, List<String>> children = new HashMap<>();
		HashSet<String> topLevels = new HashSet<>();
		while (entries.hasMoreElements()) {
			ZipEntry entry = entries.nextElement();
			long size = entry.getCompressedSize();
			Preconditions.checkState(size >= 0l);
			
			String name = entry.getName();
			if (name.endsWith("/"))
				name = name.substring(0, name.length()-1);
			Preconditions.checkState(!pathSizes.containsKey(name));
			if (name.contains("/")) {
				String[] split = name.split("/");
				System.out.println(name+" ["+split.length+"]");
				Preconditions.checkState(split.length > 1);
				String prevRunning = null;
				String running = null;
				for (int i=0; i<split.length-1; i++) {
					if (running == null)
						running = "";
					else
						running += "/";
					running += split[i];
					System.out.println("\t"+i+". "+running+"; repeat ? "+pathSizes.containsKey(running));
					if (pathSizes.containsKey(running)) {
//						System.out.println(running+" is a duplicate");
						pathSizes.put(running, pathSizes.get(running)+size);
					} else {
						// first time encoutnered
						System.out.println("Just encountered "+running);
						if (i == 0) {
							topLevels.add(running);
						} else {
							// register this as a child
							Preconditions.checkNotNull(prevRunning);
							List<String> siblings = children.get(prevRunning);
							if (siblings == null) {
								siblings = new ArrayList<>();
								children.put(prevRunning, siblings);
								System.out.println("Creating children for "+prevRunning);
							} else {
								System.out.println("Adding new child for "+prevRunning);
							}
							siblings.add(running);
						}
						pathSizes.put(running, size);
					}
					prevRunning = running;
				}
				List<String> siblings = children.get(prevRunning);
				if (siblings == null) {
					siblings = new ArrayList<>();
					children.put(prevRunning, siblings);
					System.out.println("Creating children for "+prevRunning);
				} else {
					System.out.println("Adding new child for "+prevRunning);
				}
				siblings.add(name);
			} else {
				topLevels.add(name);
			}
			pathSizes.put(name, size);
		}
		
		Map<String, Long> topSizes = new HashMap<>();
		for (String name : topLevels)
			topSizes.put(name, pathSizes.get(name));

		long expandBytes = mbToBytes(minSizeToExpand);
		long printBytes = mbToBytes(minSizeToPrint);
		
		long totalSize = 0l;
		for (String topName : ComparablePairing.getSortedData(topSizes)) {
			long bytes = pathSizes.get(topName);
			totalSize += bytes;
			System.out.println(sizeDF.format(bytesToMB(bytes))+" MB:\t"+topName);
			
			printChildrenRecursive(topName, bytes, expandBytes, printBytes, pathSizes, children, 1);
			System.out.println();
		}
		
		System.out.println("TOTAL:\t"+sizeDF.format(bytesToMB(totalSize))+" MB");
		
		zip.close();
	}
	
	private static long mbToBytes(double mb) {
		return (long)(mb*1024*1024);
	}
	
	private static double bytesToMB(long bytes) {
		return (double)bytes/(1024*1024);
	}
	
	private static final DecimalFormat sizeDF = new DecimalFormat("0.000");
	
	private static void printChildrenRecursive(String parent, long parentSize, long minExpandBytes, long minPrintBytes,
			Map<String, Long> pathSizes, Map<String, List<String>> children, int indents) {
		List<String> myChildren = children.get(parent);
		if (myChildren == null) {
//			System.out.println("No children of "+parent);
			return;
		}
//		System.out.println(myChildren.size()+" children of "+parent);
		List<Long> childSizes = new ArrayList<>(myChildren.size());
		for (String child : myChildren)
			childSizes.add(pathSizes.get(child));
		myChildren = ComparablePairing.getSortedData(childSizes, myChildren);
		Collections.reverse(myChildren);
		for (String child : myChildren) {
			long size = pathSizes.get(child);
			if (size >= minPrintBytes) {
				for (int i=0; i<indents; i++)
					System.out.print("\t");
				System.out.println(sizeDF.format(bytesToMB(size))+" MB:\t"+child);
				if (size > minExpandBytes)
					printChildrenRecursive(child, size, minExpandBytes, minPrintBytes, pathSizes, children, indents+1);
			}
		}
	}

}
