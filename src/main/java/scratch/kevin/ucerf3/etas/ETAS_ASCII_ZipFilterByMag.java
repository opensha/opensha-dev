package scratch.kevin.ucerf3.etas;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;

public class ETAS_ASCII_ZipFilterByMag {

	public static void main(String[] args) throws ZipException, IOException {
		File input = new File("/tmp/results_ascii.zip");
		File output = new File("/tmp/results_ascii_m5.zip");
		double minMag = 5d;
		
		ZipFile zip = new ZipFile(input);
		ZipOutputStream zout = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(output)));
		
		Enumeration<? extends ZipEntry> entries = zip.entries();
		int BUFFER = 8192;
		zout.setMethod(ZipOutputStream.DEFLATED);
		byte data[] = new byte[BUFFER];
		
		while (entries.hasMoreElements()) {
			ZipEntry entry = entries.nextElement();
			if (entry.getName().endsWith("simulatedEvents.txt")) {
				System.out.println("Filtering "+entry.getName());
				ETAS_Catalog cat = ETAS_CatalogIO.loadCatalog(new BufferedInputStream(zip.getInputStream(entry)), minMag);
				System.out.println("Loaded "+cat.size()+" events");
				ZipEntry newEntry = new ZipEntry(entry.getName());
				zip.getInputStream(newEntry);
//				ETAS_CatalogIO.getEventFileLine(null)
				ByteArrayOutputStream baos = new ByteArrayOutputStream();
				BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(baos));
				ETAS_CatalogIO.writeEventDataToWriter(writer, cat);
				writer.flush();
				baos.flush();
				byte[] bytes = baos.toByteArray();
				zout.putNextEntry(newEntry);
				System.out.println("writing "+bytes.length+" bytes");
				zout.write(bytes);
				zout.closeEntry();
			} else {
				System.out.println("Copying "+entry.getName());
				zout.putNextEntry(entry);
				BufferedInputStream origin = new BufferedInputStream(zip.getInputStream(entry));
				int count;
				while ( (count = origin.read(data, 0,
						BUFFER)) != -1) {
					zout.write(data, 0, count);
				}
				origin.close();
				zout.closeEntry();
			}
		}
		
		zout.close();
		
	}

}
