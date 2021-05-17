package scratch.kevin.cfm.excelParse;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.kevin.cfm.db.CFMFault;
import scratch.kevin.cfm.db.CFMFaultName;
import scratch.kevin.cfm.db.CFMFaultSection;
import scratch.kevin.cfm.db.CFMFaultStrand;
import scratch.kevin.cfm.db.CFMFaultSystem;
import scratch.kevin.cfm.db.CFMFaultZone;

public class CFM5ExcelParser {
	
	// column name/index
	// A B C D E F G H I J K  L  M  N  O  P  Q
	// 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16

	public static void main(String[] args) throws FileNotFoundException, IOException {
		File xlsFile = new File("/home/kevin/OpenSHA/cfm/CFM5-short-Fault_ID-naming.xls");
		POIFSFileSystem fs = new POIFSFileSystem(new FileInputStream(xlsFile));
		
		HSSFWorkbook wb = new HSSFWorkbook(fs);
		HSSFSheet sheet = wb.getSheetAt(0);
		
		int startRow = 3;
		int nameCol = 3; // D
		int systemStartCol = 4; // E
		int zoneStartCol = 7; // H
		int sectionStartCol = 10; // K
		int nameStartCol = 13; // L
		int strandStartCol = 16; // Q
		
		CFMFaultSystem curSystem = null;
		CFMFaultZone curZone = null;
		CFMFaultSection curSection = null;
		CFMFaultName curName = null;
		CFMFaultStrand curStrand = null;
		
		List<CFMFault> faults = Lists.newArrayList();
		
		for (int r=startRow; r<=sheet.getLastRowNum(); r++) {
			HSSFRow row = sheet.getRow(r);
			
			HSSFCell nameCell = row.getCell(nameCol);
			if (nameCell == null)
				continue;
			String objectName = nameCell.getStringCellValue();
			
			if (hasNamedGroup(row, systemStartCol)) {
				CFMFaultSystem system = new CFMFaultSystem(getGroupID(row, systemStartCol),
						getGroupName(row, systemStartCol), getGroupShortName(row, systemStartCol));
				if (curSystem == null || !system.equals(curSystem)) {
					curSystem = system;
					curZone = null;
					curSection = null;
					curName = null;
					curStrand = null;
				}
			}
			
			if (hasNamedGroup(row, zoneStartCol)) {
				CFMFaultZone zone = new CFMFaultZone(getGroupID(row, zoneStartCol),
						getGroupName(row, zoneStartCol), getGroupShortName(row, zoneStartCol));
				if (curZone == null || !zone.equals(curZone)) {
					curZone = zone;
					curSection = null;
					curName = null;
					curStrand = null;
				}
			}
			
			if (hasNamedGroup(row, sectionStartCol)) {
				CFMFaultSection section = new CFMFaultSection(getGroupID(row, sectionStartCol),
						getGroupName(row, sectionStartCol), getGroupShortName(row, sectionStartCol));
				if (curSection == null || !section.equals(curSection)) {
					curSection = section;
					curName = null;
					curStrand = null;
				}
			}
			
			if (hasNamedGroup(row, nameStartCol)) {
				CFMFaultName name = new CFMFaultName(getGroupID(row, nameStartCol),
						getGroupName(row, nameStartCol), getGroupShortName(row, nameStartCol));
				if (name == null || !name.equals(curName)) {
					curName = name;
					curStrand = null;
				}
			}
			
			if (hasNamedGroup(row, strandStartCol)) {
//				CFMFaultStrand strand = new CFMFaultStrand(getGroupID(row, strandStartCol),
//						getGroupName(row, strandStartCol), getGroupShortName(row, strandStartCol));
//				if (curStrand == null || !strand.equals(curStrand))
//					curStrand = strand;
			}
			
			Preconditions.checkNotNull(curSystem);
			Preconditions.checkNotNull(curZone);
			Preconditions.checkNotNull(curSection);
			Preconditions.checkNotNull(curName);
			String autoPrefix = curSystem.getShortName()+"-"+curZone.getShortName()+"-"+curSection.getShortName();
			
			if (!objectName.startsWith(autoPrefix)) {
				System.out.println("\t***SKIPPING: "+objectName);
				continue;
			}
			
			CFMFault fault = new CFMFault(objectName, curSystem, curZone, curSection, curName, curStrand);
			System.out.println(fault);
			
			faults.add(fault);
		}
	}
	
	private static boolean hasNamedGroup(HSSFRow row, int groupStartCol) {
		HSSFCell cell = row.getCell(groupStartCol);
		return cell != null && !getGroupName(row, groupStartCol).isEmpty();
//		try {
//			return cell != null && getGroupID(row, groupStartCol) >= 0;
//		} catch (NumberFormatException e) {
//			return false;
//		}
	}
	
	private static int getGroupID(HSSFRow row, int groupStartCol) {
		return -1;
//		HSSFCell cell = row.getCell(groupStartCol);
//		Preconditions.checkNotNull(cell);
//		if (cell.getCellType() == HSSFCell.CELL_TYPE_NUMERIC)
//			return (int)cell.getNumericCellValue();
//		else
//			return Integer.parseInt(cell.getStringCellValue());
		
		
//		Preconditions.checkState(cell != null && cell.getCellType() == HSSFCell.CELL_TYPE_NUMERIC);
//		return (int)cell.getNumericCellValue();
	}
	
	private static String getGroupName(HSSFRow row, int groupStartCol) {
		HSSFCell cell = row.getCell(groupStartCol+1);
		if (cell == null || cell.getCellType() == HSSFCell.CELL_TYPE_BLANK)
			return "";  // blank ok for name
		Preconditions.checkState(cell.getCellType() == HSSFCell.CELL_TYPE_STRING);
		return cell.getStringCellValue();
	}
	
	private static String getGroupShortName(HSSFRow row, int groupStartCol) {
		HSSFCell cell = row.getCell(groupStartCol+2);
		if (cell == null || cell.getCellType() == HSSFCell.CELL_TYPE_BLANK)
			return getGroupName(row, groupStartCol).replaceAll(" ", "");
		Preconditions.checkNotNull(cell);
		if (cell.getCellType() == HSSFCell.CELL_TYPE_NUMERIC) {
			double val = cell.getNumericCellValue();
			if (val == Math.round(val))
				return (int)val+"";
			else
				return val+"";
		}
		Preconditions.checkState(cell.getCellType() == HSSFCell.CELL_TYPE_STRING);
		if (cell.getStringCellValue().isEmpty())
			return getGroupName(row, groupStartCol).replaceAll(" ", "");
		return cell.getStringCellValue();
	}

}
