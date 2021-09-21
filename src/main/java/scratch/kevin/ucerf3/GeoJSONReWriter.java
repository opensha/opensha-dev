package scratch.kevin.ucerf3;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.Geometry.DepthSerializationType;

import com.google.common.io.Files;
import com.google.gson.Gson;

public class GeoJSONReWriter {

	public static void main(String[] args) throws ZipException, IOException {
		File mainDir = new File("/home/kevin/markdown/inversions");
		for (File runDir : mainDir.listFiles()) {
			if (!runDir.isDirectory())
				continue;
			System.out.println(runDir.getName());
			for (File file : new File[] {new File(runDir, "rupture_set.zip"), new File(runDir, "solution.zip")}) {
				if (file.exists()) {
					System.out.println("\tProcessing "+file.getName());
					ZipFile zip = new ZipFile(file);
					
					File outTmpFile = new File(file.getAbsolutePath()+".tmp");
					
					ZipOutputStream zout = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(outTmpFile)));
					
					Enumeration<? extends ZipEntry> entries = zip.entries();
					while (entries.hasMoreElements()) {
						ZipEntry entry = entries.nextElement();
						
						boolean copyOver = true;
						
						if (entry.getName().endsWith("fault_sections.geojson")) {
							// it's  a match, modify
							System.out.println("\tConverting "+entry.getName());
							
//							FeatureCollection.
							Gson inGSON = FeatureCollection.buildGson(DepthSerializationType.ELEVATION_M);
							Gson outGSON = FeatureCollection.buildGson(DepthSerializationType.DEPTH_KM);
							
							BufferedReader read = new BufferedReader(new InputStreamReader(zip.getInputStream(entry)));
							
							try {
								FeatureCollection features = inGSON.fromJson(read, FeatureCollection.class);
								read.close();
								
								zout.flush();
								zout.putNextEntry(new ZipEntry(entry.getName()));
								BufferedWriter write = new BufferedWriter(new OutputStreamWriter(zout));
								outGSON.toJson(features, FeatureCollection.class, write);
								write.flush();
								zout.closeEntry();
								copyOver = false;
							} catch (Exception e) {
								e.printStackTrace();
							}
						}
						
						if (copyOver) {
							// copy it over
							zout.putNextEntry(new ZipEntry(entry.getName()));
							
							BufferedInputStream bin = new BufferedInputStream(zip.getInputStream(entry));
							bin.transferTo(zout);
							zout.flush();
							
							zout.closeEntry();
						}
					}
					zout.close();
					zip.close();
					System.out.println("\tDone. moving "+outTmpFile.getName()+" to "+file.getName());
					Files.move(file, new File(file.getAbsoluteFile()+".orig"));
					Files.move(outTmpFile, file);
				}
			}
		}
	}

}
