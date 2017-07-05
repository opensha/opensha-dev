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


/**
 * This class provides a simple example of generating and sending a product.
 * A product is a bundle of information including properties (key-value pairs)
 * as well as contents (files/blobs). Read more about products here:
 *
 * http://ehppdl1.cr.usgs.gov/userguide/products/index.html
 *
 * This example is generalized for how to send products programmatically and
 * does not reflect the expected properties or contents that might need to be
 * added to an actual OAF product.
 *
 * ```
 * Final OAF product properties, contents and structure requires discussion.
 * ```
 *
 * For OAF, there are currently three types of products that the event pages
 * will look for:
 *
 * 1) type = "oaf"
 *      This is the primary OAF product and should include forecast
 *      information. This product is used to generate the forecast table and
 *      fill in the boilerplate text on the OAF section of the event page. This
 *      product must be sent in order to generate an OAF section on the event
 *      page.
 *
 * 2) type = "oaf-header"
 *      This product allows free-form introductory text to be added to the
 *      top of the OAF section of the event page (above the tabs). No special
 *      styling is applied to this product when rendered, but the product
 *      content (inline only) may be set to use a type of "text/html" to allow
 *      for some simple styling.
 *      See: http://ehppdl1.cr.usgs.gov/userguide/products/text.html
 *
 * 3) type = "oaf-link"
 *      This product allows for links to be added to the OAF section of the
 *      event page. A link is probably added in order to direct visitors to
 *      additional information, like local emergency reponse organizations.
 *      See: http://ehppdl1.cr.usgs.gov/userguide/products/link.html
 *
 * We can discuss if additional product types are necessary, as well as discuss
 * the final form for the base "oaf" product type in order to ensure the
 * OAF section on the event page can be generated as desired.
 *
 * Contact
 *  - Email: "Eric Martinez" <emartinez@usgs.gov>
 *  - Phone: 303.273.8653
 */
public class OafSendExample {
	/**
	 * Simple example to send a product.
	 *
	 * Step 1: Create a product {gov.usgs.earthquake.product.Product}
	 * Step 2: Create a sender  {gov.usgs.earthquake.distribution.SocketProductSender}
	 * Step 3: Use sender to send the product
	 *
	 */
	public static void main (String [] args) throws Exception {
		if ("".isEmpty())
			throw new RuntimeException("Probably shouldn't actually run this example from Eric");
		Product product;
		ProductSender sender;

		product = createProduct();
		sender = createSender();

		// For production, you will need to sign the product
		// File privateKey = new File("~/.ssh/id_rsa"); // OpenSSH private key file
		// product.sign(CryptoUtils.readOpenSSHPrivateKey(StreamUtils.readStream(
		//     StreamUtils.getInputStream(privateKey)), null));

		// ... may throw exceptions to handle ...
		sender.sendProduct(product);
	}


	/**
	 * If you have data structures/objects in runtime memory and want to attach
	 * them to the product, this is a good way to do that.
	 *
	 * @param product
	 *     The product on to which the content should be attached.
	 */
	public static void attachByteContentToProduct (Product product) {
		ByteContent diskContent, inlineContent;
		Map<String, Content> contents;

		contents = product.getContents();

		// This content will be stored on disk in ComCat. For the event pages
		// to access this content, they must make a separate HTTP request (which
		// can slow page rendering). The only thing that makes this "disk"
		// content is that it receives a non-empty content-path when
		// it is added to the contents map. It has nothing to do with the
		// origination of the content itself... It is perfectly fine for this
		// content to be binary blobs if necessary.
		diskContent = new ByteContent("Here is some disk content".getBytes());
		diskContent.setContentType(/*String*/ null); // `null` implies "text/plain"
		diskContent.setLastModified(/*Date*/ null); // `null` implies NOW
		contents.put("disk-content.txt", diskContent);

		// This content will be embedded into the event details feed from ComCat
		// and is immediately available to the event page. The only thing that
		// makes this "inline" content is that the content-path is an empty string.
		// Content added as inline should be text-based and _not_ binary blobs.
		inlineContent = new ByteContent(
				"{\"key\": \"Here is some inline content\"}".getBytes());
		inlineContent.setContentType("application/json");
		inlineContent.setLastModified(/*Date*/ null); // `null` implies now
		contents.put("", inlineContent); // <-- Inline content because empty string

		// In general, you can serialize any Java object to a String and then
		// use String.getBytes() to generate ByteContent for a product.
	}

	/**
	 * Attach contents to the product. Internally a product stores contents as
	 * a map from content-path to Content object. The content-path indicates
	 * the relative path/filename of the content. Note, this applies to both
	 * byte-based (runtime objects) and file-based Content.
	 *
	 * You may create _one_ simple "inline" Content object by specifying the
	 * empty string "" as the content-path. This inline content is more readily
	 * accessibly to the event pages as it avoids an additional HTTP request.
	 *
	 * @param product
	 *     The product on to which the content should be attached.
	 */
	public static void attachContentToProduct (Product product)
			throws Exception {
		attachByteContentToProduct(product);
		attachFileContentToProduct(product);
	}

	/**
	 * Attaches files stored locally on disk to the product contents.
	 *
	 * @param product
	 *     The product on to which the content should be attached.
	 */
	public static void attachFileContentToProduct (Product product)
			throws Exception {
		FileContent fileContent;
		Map<String, Content> contents;

		contents = product.getContents();

		// Add a single file to the contents
		fileContent = new FileContent(new File("./contents.xml"));
		contents.put("contents.xml", fileContent);

		// Add all files in a directory to the contents. Note, this may
		// throw an exception...
		contents.putAll(FileContent.getDirectoryContents(
				new File("./test-content")));
	}

	/**
	 * Attaches properties (key-value pairs) to a product.
	 *
	 * @param product
	 *     The product on to which the properties should be attached.
	 */
	public static void attachPropertiesToProduct (Product product) {
		Map<String, String> properties;

		// Properties is a reference to the product's internal properties hash,
		// so updates made here are immediately reflected on the product properties
		properties = product.getProperties();


		// Typically, you always want to set the eventsource and eventsourcecode
		// properties. These are the strongest way to ensure the product being
		// sent will associate to the correct event once it reaches ComCat

		// Note: We could also use the NC information "nc" and "72689331"
		// or any other eventid that references the same event that actually
		// occurred...
		properties.put(Product.EVENTSOURCE_PROPERTY, "us");
		properties.put(Product.EVENTSOURCECODE_PROPERTY, "10006jv5");


		// You may also want to set these other event-related properties
		// which are used as a fall-back for association purposes. There are pros
		// and cons to doing so that I won't discuss here. Dealer's choice!

		// properties.put(Product.EVENTTIME_PROPERTY, "2016-09-03T03:27:57.170UTC");
		// properties.put(Product.MAGNITUDE_PROPERTY, "5.7");
		// properties.put(Product.LATITUDE_PROPERTY, "40.406");
		// properties.put(Product.LONGITUDE_PROPERTY, "-125.469");

		// properties.put(Product.DEPTH_PROPERTY, "5.6");

		// Note: Product.VERSION_PROPERTY does not indicate a "version" within
		// ComCat, but is really just an identifier for senders/network operators.
		// Typically it is not recommended to send this property.

		// properties.put(Proudct.VERSION_PROPERTY, "01");


		// This property is interpreted on the event page and dicates whether
		// an orange "X" or green checkmark is displayed at the top of the
		// product page. This status icon informs the user if this is an
		// automatically-generated forecast, or if a scientist manually generated
		// it. If this property is omitted, it is assumed the product was
		// automatically generated.
		properties.put("review-status", "reviewed"); // "reviewed" or "automatic"


		// Also add properties specific to your product-type, in this case,
		// things that are specific to the aftershock forecast.
		// We may want to discuss what specifically you add here...

		properties.put("wcradius", "167.55");
		// properties.put("key", "value"); // ... as needed/appropriate ...
	}

	/**
	 *
	 * @return
	 *     The product that was created
	 */
	public static Product createProduct () throws Exception {
		Product product;
		ProductId productId;

		// Need an ID for the product. The ID is for the product only and does
		// not necessarily reflect "event" properties, however, in practice
		// they are typically quite similar

		// A product ID uniquely identifies a product based on "source", "type",
		// and "code" properties. An "updateTime" property indicates a specific
		// version of a product based on when the product was created. By default,
		// newer versions of the same product (source/type/code) will replace
		// older versions.

		// Read more about source, type, code here:
		// http://ehppdl1.cr.usgs.gov/userguide/products/index.html
		productId = new ProductId(
				"some-source",   // source
				"some-type",     // type
				"some-code"      // code
				// updateTime - Defaults to NOW if not specificied, typically, do not
				//              specify this as a best practice
				);

		// Use Product.STATUS_UPDATE to create or update a product.
		// Use Product.STATUS_DELETE to remove this product.
		product = new Product(
				productId,
				Product.STATUS_UPDATE
				);

		// This is a legacy vestige that does not get used anymore, but we need
		// to add it otherwise sending will fail. Anything will do ...
		product.setTrackerURL(new URL("http://www.google.com/"));


		attachPropertiesToProduct(product); // simply key-value properties
		attachContentToProduct(product);    // file- or byte-based content

		return product;
	}

	/**
	 * Creates a sender to push products into PDL.
	 *
	 */
	public static ProductSender createSender () {
		SocketProductSender sender;

		// SocketProductSenders send directly to a PDL HUB and do not introduce
		// any polling latency.
		sender = new SocketProductSender("pdldevel.cr.usgs.gov", 11235);

		// If product consists primarily of binary data, set this option `true`
		// to accelerate distribution.
		sender.setBinaryFormat(false);

		// If product consists primarily of text data, set this option `true`
		// to accelerate distribution.
		sender.setEnableDeflate(true);

		// ^^ Note ^^ Typically do not set both of the above options to `true` as
		//            binary content doesn't compress efficiently but adds
		//            processing overhead.


		return sender;
	}
}