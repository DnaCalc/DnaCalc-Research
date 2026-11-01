// Original research adapter; the called Apache POI sources retain their own licence.
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import org.apache.poi.ss.util.NumberToTextConverter;

public final class NumberTextProbe {
    public static void main(String[] args) throws Exception {
        try (BufferedReader input = new BufferedReader(new InputStreamReader(System.in, StandardCharsets.US_ASCII))) {
            String line;
            while ((line = input.readLine()) != null) {
                long bits = Long.parseUnsignedLong(line, 16);
                System.out.println(line + "\t" + NumberToTextConverter.toText(Double.longBitsToDouble(bits)));
            }
        }
    }
}
