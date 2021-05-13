package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;

public class CommonsImportCheck {
	
	public static Map<String, List<String>> getImportsOutside(File fromPackageDir,
			String matchesPrefix, String... whitelistPrefixes) throws IOException {
		Map<String, List<String>> imports = Maps.newHashMap();
		for (File sourceFile : fromPackageDir.listFiles()) {
			if (sourceFile.isDirectory()) {
				imports.putAll(getImportsOutside(sourceFile, matchesPrefix, whitelistPrefixes));
				continue;
			}
			String fName = sourceFile.getName();
			if (!fName.toLowerCase().endsWith(".java"))
				continue;
			List<String> myImports = null;
			String myName = fName.substring(0, fName.toLowerCase().indexOf(".java"));
			String myPackage = null;
			
			for (String line : Files.readLines(sourceFile, Charset.defaultCharset())) {
				if (myPackage == null && line.trim().startsWith("package")) {
					// parse the package
					myPackage = line.substring(line.indexOf("package")+"package".length());
					myPackage = myPackage.substring(0, myPackage.indexOf(";")).trim();
				}
				if (line.trim().startsWith("import")) {
					String myImport = line.substring(line.indexOf("import")+"import".length());
					myImport = myImport.substring(0, myImport.indexOf(";")).trim();
//					System.out.println("import line: "+line+" parsed to "+myImport);
					if (myImport.startsWith(matchesPrefix)) {
						// potential bad import
						boolean whitelist = false;
						for (String whitelistPrefix : whitelistPrefixes) {
							if (myImport.startsWith(whitelistPrefix)) {
								whitelist = true;
								break;
							}
						}
						if (!whitelist) {
							// bad import
							if (myImports == null)
								myImports = Lists.newArrayList();
							myImports.add(myImport);
						}
					}
				}
			}
			if (myImports != null)
				imports.put(myPackage+"."+myName, myImports);
		}
		
		return imports;
	}

	public static void main(String[] args) throws IOException {
		File commonsDir = new File("src"+File.separator+"org"+File.separator+"opensha"+File.separator+"commons");
		Map<String, List<String>> imports = getImportsOutside(commonsDir, "org.opensha", "org.opensha.commons");
		for (String className : imports.keySet()) {
			System.out.println(className);
			for (String imported : imports.get(className))
				System.out.println("\t"+imported);
		}
		System.out.println(imports.size()+" offenders");
	}

}
