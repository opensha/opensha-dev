package scratch.aftershockStatistics.pdl;

import gov.usgs.earthquake.distribution.ProductSender;
import gov.usgs.earthquake.distribution.SocketProductSender;

import gov.usgs.earthquake.product.ByteContent;
import gov.usgs.earthquake.product.Content;
import gov.usgs.earthquake.product.FileContent;
import gov.usgs.earthquake.product.Product;
import gov.usgs.earthquake.product.ProductId;

import gov.usgs.util.CryptoUtils;
import gov.usgs.util.StreamUtils;

import java.io.File;
import java.net.URL;
import java.util.Map;
import java.util.List;



/**
 * Code to construct the PDL product of type "oaf".
 * Author: Michael Barall 06/06/2018.
 */
public class PDLProductBuilderOaf {




	/**
	 * Create the PDL product.
	 * @param eventID = Event ID, used as the "code" for the product (for example, "us10006jv5").
	 * @param eventNetwork = Network identifier for the event (for example, "us").
	 * @param eventCode = Network code for the event (for example, "10006jv5").
	 * @param isReviewed = True if this product has been reviewed.
	 * @param jsonText = JSON text that contains the product.
	 * @param modifiedTime = Modification time, in milliseconds since the epoch, or 0L if none.
	 * Note: At present, modifiedTime is ignored, and the time is always set to "now".
	 * @return
	 *     The product that was created
	 */
	public static Product createProduct (String eventID, String eventNetwork, String eventCode, boolean isReviewed, String jsonText, long modifiedTime) throws Exception {
		Product product;
		ProductId productId;

//		// Need an ID for the product. The ID is for the product only and does
//		// not necessarily reflect "event" properties, however, in practice
//		// they are typically quite similar
//
//		// A product ID uniquely identifies a product based on "source", "type",
//		// and "code" properties. An "updateTime" property indicates a specific
//		// version of a product based on when the product was created. By default,
//		// newer versions of the same product (source/type/code) will replace
//		// older versions.
//
//		// Read more about source, type, code here:
//		// http://ehppdl1.cr.usgs.gov/userguide/products/index.html
//		productId = new ProductId(
//				"some-source",   // source
//				"some-type",     // type
//				"some-code"      // code
//				// updateTime - Defaults to NOW if not specificied, typically, do not
//				//              specify this as a best practice
//				);
//
//		// Use Product.STATUS_UPDATE to create or update a product.
//		// Use Product.STATUS_DELETE to remove this product.
//		product = new Product(
//				productId,
//				Product.STATUS_UPDATE
//				);
//
//		// This is a legacy vestige that does not get used anymore, but we need
//		// to add it otherwise sending will fail. Anything will do ...
//		product.setTrackerURL(new URL("http://www.google.com/"));
//
//
//		attachPropertiesToProduct(product); // simply key-value properties
//		attachContentToProduct(product);    // file- or byte-based content


		// Check arguments

		if (eventID == null || eventID.isEmpty()) {
			throw new IllegalArgumentException ("PDLProductBuilderOaf: Event ID is not specified");
		}
		if (eventNetwork == null || eventNetwork.isEmpty()) {
			throw new IllegalArgumentException ("PDLProductBuilderOaf: Event network identifier is not specified");
		}
		if (eventCode == null || eventCode.isEmpty()) {
			throw new IllegalArgumentException ("PDLProductBuilderOaf: Event network code is not specified");
		}
		if (jsonText == null || jsonText.isEmpty()) {
			throw new IllegalArgumentException ("PDLProductBuilderOaf: JSON text for product content is not specified");
		}

		// Announce it

		System.out.println ("Constructing PDL oaf product: event ID = " + eventID + ", network = " + eventNetwork + ", code = " + eventCode);

		// Our PDL source, type, and code

		String source = "us";
		String type = "oaf";
		String code = eventID;

		// Construct the product ID

		productId = new ProductId(
				source,   // source
				type,     // type
				code      // code
				// updateTime - Defaults to NOW if not specificied, typically, do not
				//              specify this as a best practice
				);

		// Construct the product

		product = new Product(
				productId,
				Product.STATUS_UPDATE
				);

		// Legacy requirement

		product.setTrackerURL(new URL("http://www.google.com/"));

		// Attach properties

		attachPropertiesToProduct (product, eventNetwork, eventCode, isReviewed);

		// Attach content

		attachByteContentToProduct (product, jsonText, modifiedTime);

		// Return the constructed product

		return product;
	}




	/**
	 * Create a PDL deletion product.
	 * @param eventID = Event ID, used as the "code" for the product (for example, "us10006jv5").
	 * @param eventNetwork = Network identifier for the event (for example, "us").
	 * @param eventCode = Network code for the event (for example, "10006jv5").
	 * @param isReviewed = True if this product has been reviewed.
	 * @param modifiedTime = Modification time, in milliseconds since the epoch, or 0L if none.
	 * Note: At present, modifiedTime is ignored, and the time is always set to "now".
	 * @return
	 *     The product that was created
	 */
	public static Product createDeletionProduct (String eventID, String eventNetwork, String eventCode, boolean isReviewed, long modifiedTime) throws Exception {
		Product product;
		ProductId productId;

//		// Need an ID for the product. The ID is for the product only and does
//		// not necessarily reflect "event" properties, however, in practice
//		// they are typically quite similar
//
//		// A product ID uniquely identifies a product based on "source", "type",
//		// and "code" properties. An "updateTime" property indicates a specific
//		// version of a product based on when the product was created. By default,
//		// newer versions of the same product (source/type/code) will replace
//		// older versions.
//
//		// Read more about source, type, code here:
//		// http://ehppdl1.cr.usgs.gov/userguide/products/index.html
//		productId = new ProductId(
//				"some-source",   // source
//				"some-type",     // type
//				"some-code"      // code
//				// updateTime - Defaults to NOW if not specificied, typically, do not
//				//              specify this as a best practice
//				);
//
//		// Use Product.STATUS_UPDATE to create or update a product.
//		// Use Product.STATUS_DELETE to remove this product.
//		product = new Product(
//				productId,
//				Product.STATUS_UPDATE
//				);
//
//		// This is a legacy vestige that does not get used anymore, but we need
//		// to add it otherwise sending will fail. Anything will do ...
//		product.setTrackerURL(new URL("http://www.google.com/"));
//
//
//		attachPropertiesToProduct(product); // simply key-value properties
//		attachContentToProduct(product);    // file- or byte-based content


		// Check arguments

		if (eventID == null || eventID.isEmpty()) {
			throw new IllegalArgumentException ("PDLProductBuilderOaf: Event ID is not specified");
		}
		if (eventNetwork == null || eventNetwork.isEmpty()) {
			throw new IllegalArgumentException ("PDLProductBuilderOaf: Event network identifier is not specified");
		}
		if (eventCode == null || eventCode.isEmpty()) {
			throw new IllegalArgumentException ("PDLProductBuilderOaf: Event network code is not specified");
		}

		// Announce it

		System.out.println ("Constructing PDL deletion oaf product: event ID = " + eventID + ", network = " + eventNetwork + ", code = " + eventCode);

		// Our PDL source, type, and code

		String source = "us";
		String type = "oaf";
		String code = eventID;

		// Construct the product ID

		productId = new ProductId(
				source,   // source
				type,     // type
				code      // code
				// updateTime - Defaults to NOW if not specificied, typically, do not
				//              specify this as a best practice
				);

		// Construct the product

		product = new Product(
				productId,
				Product.STATUS_DELETE
				);

		// Legacy requirement

		product.setTrackerURL(new URL("http://www.google.com/"));

		// Attach properties

		attachPropertiesToProduct (product, eventNetwork, eventCode, isReviewed);

		// Return the constructed product

		return product;
	}




	/**
	 * If you have data structures/objects in runtime memory and want to attach
	 * them to the product, this is a good way to do that.
	 *
	 * @param product
	 *     The product on to which the content should be attached.
	 * @param jsonText = JSON text that contains the product.
	 * @param modifiedTime = Modification time, in milliseconds since the epoch, or 0L if none.
	 * Note: At present, modifiedTime is ignored, and the time is always set to "now".
	 */
	public static void attachByteContentToProduct (Product product, String jsonText, long modifiedTime) {
		Map<String, Content> contents;

		contents = product.getContents();

//		// This content will be stored on disk in ComCat. For the event pages
//		// to access this content, they must make a separate HTTP request (which
//		// can slow page rendering). The only thing that makes this "disk"
//		// content is that it receives a non-empty content-path when
//		// it is added to the contents map. It has nothing to do with the
//		// origination of the content itself... It is perfectly fine for this
//		// content to be binary blobs if necessary.
//		ByteContent diskContent = new ByteContent("Here is some disk content".getBytes());
//		diskContent.setContentType(/*String*/ null); // `null` implies "text/plain"
//		diskContent.setLastModified(/*Date*/ null); // `null` implies NOW
//		contents.put("disk-content.txt", diskContent);
//
//		// This content will be embedded into the event details feed from ComCat
//		// and is immediately available to the event page. The only thing that
//		// makes this "inline" content is that the content-path is an empty string.
//		// Content added as inline should be text-based and _not_ binary blobs.
//		ByteContent inlineContent = new ByteContent(
//				"{\"key\": \"Here is some inline content\"}".getBytes());
//		inlineContent.setContentType("application/json");
//		inlineContent.setLastModified(/*Date*/ null); // `null` implies now
//		contents.put("", inlineContent); // <-- Inline content because empty string
//
//		// In general, you can serialize any Java object to a String and then
//		// use String.getBytes() to generate ByteContent for a product.


		// Check arguments

		if (jsonText == null || jsonText.isEmpty()) {
			throw new IllegalArgumentException ("PDLProductBuilderOaf: JSON text for product content is not specified");
		}

		// Attach our JSON text as inline content

		ByteContent inlineContent = new ByteContent(jsonText.getBytes());
		inlineContent.setContentType("application/json");
		inlineContent.setLastModified(/*Date*/ null); // `null` implies now
		contents.put("", inlineContent); // <-- Inline content because empty string

		return;
	}




	/**
	 * Attaches files stored locally on disk to the product contents.
	 *
	 * @param product
	 *     The product on to which the content should be attached.
	 */
	public static void attachFileContentToProduct (Product product) throws Exception {
		Map<String, Content> contents;

		contents = product.getContents();

//		// Add a single file to the contents
//		FileContent fileContent = new FileContent(new File("./contents.xml"));
//		contents.put("contents.xml", fileContent);
//
//		// Add all files in a directory to the contents. Note, this may
//		// throw an exception...
//		contents.putAll(FileContent.getDirectoryContents(
//				new File("./test-content")));

		return;
	}




	/**
	 * Attaches properties (key-value pairs) to a product.
	 *
	 * @param product
	 *     The product on to which the properties should be attached.
	 * @param eventNetwork = Network identifier for the event (for example, "us").
	 * @param eventCode = Network code for the event (for example, "10006jv5").
	 * @param isReviewed = True if this product has been reviewed.
	 */
	public static void attachPropertiesToProduct (Product product, String eventNetwork, String eventCode, boolean isReviewed) {
		Map<String, String> properties;

		// Properties is a reference to the product's internal properties hash,
		// so updates made here are immediately reflected on the product properties
		properties = product.getProperties();


//		// Typically, you always want to set the eventsource and eventsourcecode
//		// properties. These are the strongest way to ensure the product being
//		// sent will associate to the correct event once it reaches ComCat
//
//		// Note: We could also use the NC information "nc" and "72689331"
//		// or any other eventid that references the same event that actually
//		// occurred...
//		properties.put(Product.EVENTSOURCE_PROPERTY, "us");
//		properties.put(Product.EVENTSOURCECODE_PROPERTY, "10006jv5");
//
//
//		// You may also want to set these other event-related properties
//		// which are used as a fall-back for association purposes. There are pros
//		// and cons to doing so that I won't discuss here. Dealer's choice!
//
//		// properties.put(Product.EVENTTIME_PROPERTY, "2016-09-03T03:27:57.170UTC");
//		// properties.put(Product.MAGNITUDE_PROPERTY, "5.7");
//		// properties.put(Product.LATITUDE_PROPERTY, "40.406");
//		// properties.put(Product.LONGITUDE_PROPERTY, "-125.469");
//
//		// properties.put(Product.DEPTH_PROPERTY, "5.6");
//
//		// Note: Product.VERSION_PROPERTY does not indicate a "version" within
//		// ComCat, but is really just an identifier for senders/network operators.
//		// Typically it is not recommended to send this property.
//
//		// properties.put(Proudct.VERSION_PROPERTY, "01");
//
//
//		// This property is interpreted on the event page and dicates whether
//		// an orange "X" or green checkmark is displayed at the top of the
//		// product page. This status icon informs the user if this is an
//		// automatically-generated forecast, or if a scientist manually generated
//		// it. If this property is omitted, it is assumed the product was
//		// automatically generated.
//		properties.put("review-status", "reviewed"); // "reviewed" or "automatic"
//
//
//		// Also add properties specific to your product-type, in this case,
//		// things that are specific to the aftershock forecast.
//		// We may want to discuss what specifically you add here...
//
//		properties.put("wcradius", "167.55");
//		// properties.put("key", "value"); // ... as needed/appropriate ...


		// Check arguments

		if (eventNetwork == null || eventNetwork.isEmpty()) {
			throw new IllegalArgumentException ("PDLProductBuilderOaf: Event network identifier is not specified");
		}
		if (eventCode == null || eventCode.isEmpty()) {
			throw new IllegalArgumentException ("PDLProductBuilderOaf: Event network code is not specified");
		}

		// This identifies the event

		properties.put(Product.EVENTSOURCE_PROPERTY, eventNetwork);
		properties.put(Product.EVENTSOURCECODE_PROPERTY, eventCode);

		// Indicate if the event has been reviewed

		if (isReviewed) {
			properties.put("review-status", "reviewed"); // "reviewed" or "automatic"
		}

		return;
	}

}
