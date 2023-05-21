import java.io.File;
import java.io.IOException;
import java.util.Locale;
import java.util.Scanner;

/**
 *  This class provides methods for reading ints and doubles for my selection
 *  and reads in data of various types from standard input, files.
 *  <p>
 *  For additional documentation, see
 *  https://github.com/ZjuWxy/SelectionForDBPs
 *
 *  @author Wang Xinyang
 */
public final class In {
    
    private Scanner scanner;

    // assume Unicode UTF-8 encoding
    private static final String CHARSET_NAME = "UTF-8";

    
    /**
     * Read and return the next int.
     */
    public int readInt() {
        return scanner.nextInt();
    }

    /**
     * Read and return the next double.
     */
    public double readDouble() {
        return scanner.nextDouble();
    }

    /**
     * Create an input stream from a filename.
     */
    public In(String s) {
        try {
            // first try to read file from local file system
            File file = new File(s);
            if (file.exists()) {
                scanner = new Scanner(file, CHARSET_NAME);
                scanner.useLocale(LOCALE);
                return;
            }
        }
        catch (IOException ioe) {
            System.err.println(s);
        }
    }

}
