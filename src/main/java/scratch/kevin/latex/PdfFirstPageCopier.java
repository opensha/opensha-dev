package scratch.kevin.latex;

import com.google.common.io.Files;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.pdf.PdfCopy;
import com.itextpdf.text.pdf.PdfReader;
import com.itextpdf.text.pdf.PdfStamper;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashMap;

import org.apache.commons.lang3.exception.ExceptionUtils;

public class PdfFirstPageCopier {

	/**
	 * Copies only the first page of the source PDF to the destination.
	 * If the source has exactly 1 page, the full file is copied directly.
	 *
	 * Document metadata (info dictionary) is preserved.
	 */
	public static void copyFirstPage(File src, File dest) throws IOException {
		PdfReader reader = new PdfReader(src.getAbsolutePath());
		try {
			int numPages = reader.getNumberOfPages();

			if (numPages <= 1) {
				// Only one page → copy whole file
				Files.copy(src, dest);
				return;
			}

			// Extract metadata
			HashMap<String, String> info = reader.getInfo();

			// --- Step 1: Create temp file for first-page-only PDF ---
			File tmp = File.createTempFile("first_page_", ".pdf");
			tmp.deleteOnExit();

			Document document = new Document(reader.getPageSizeWithRotation(1));
			PdfCopy copy = new PdfCopy(document, new FileOutputStream(tmp));

			document.open();
			copy.addPage(copy.getImportedPage(reader, 1));
			document.close();
			copy.close();

			reader.close(); // close original reader early

			// --- Step 2: Apply metadata via PdfStamper ---
			PdfReader tmpReader = new PdfReader(tmp.getAbsolutePath());
			PdfStamper stamper = new PdfStamper(tmpReader, new FileOutputStream(dest));
			stamper.setMoreInfo(info);
			stamper.close();
			tmpReader.close();
			
			tmp.delete();

		} catch (DocumentException e) {
			throw ExceptionUtils.asRuntimeException(e);
		} finally {
			reader.close();
		}
	}
}
