import java.io.FileNotFoundException;
import java.io.PrintStream;
import static java.lang.Math.abs;

public class origin {
    // Standards
    private static int len = 6750;

    private static double limMZ = 0.01;
    private static double limTIC = 0.5;
    private static double limRT = 0.15;

    // read in data
    public static entry[] readFromTxt (String st){
        In in = new In(st);
        entry[] res = new entry[len+1];
        for (int i = 1; i <= len; i++){
            res[i] = new entry(in.readInt(),in.readDouble(),in.readDouble(),in.readDouble(),in.readDouble(),
                    in.readInt(), in.readInt());
        }
        return res;
    }

    public static void main(String[] args) throws FileNotFoundException {
        entry[] entries = readFromTxt("data.txt");
        PrintStream ps1 = new PrintStream("res.txt");
        System.setOut(ps1);

        // lLim & rLim
        double lMZ = 2.00-limMZ;
        double rMZ = 2.00+limMZ;
        double lTIC = 3.00-limTIC;
        double rTIC = 3.00+limTIC;
        System.out.println("筛选标准_1: 质荷比:  M : (M+" + lMZ + ") ~ (M+" + rMZ +") ");
        System.out.println("筛选标准_2: 丰度比:  1.0 : " + lTIC + " ~ " + rTIC);
        System.out.println("筛选标准_3: 保留时间差:  |RT_1 - RT_2| < " + limRT);
        System.out.println();
        System.out.println("可能的含氯样本: ");

        // #sample
        int counts = 1;
        // two iteration
        for (int i = 1; i < len; i++){
            double t1, t2, t3;
            for(int k = 1; k < i; k++){
                t1 = abs(entries[i].mz - entries[k].mz - 2);
                t2 = abs(entries[k].TIC / entries[i].TIC - 3);
                t3 = abs(entries[k].RT - entries[i].RT);
                // standards to meet
                if (t1 <= limMZ && t2 <= limTIC && t3 <= limRT){
                    System.out.println("样本号: " + counts++);
                    System.out.println("   索引号: " + "#" + entries[k].index + " 和 " + "#" +  entries[i].index);
                    double s1 = t1 + 2;
                    double s2 = t2 + 100/97.8;
                    System.out.println("   质荷比差: " + String.format("%.4f", s1) +
                            "   丰度比: 100.0% , " + String.format("%.1f", entries[i].TIC/entries[k].TIC * 100.0) + "%" +
                            "   保留时间差: " + String.format("%.4f", t3));
                    System.out.println("   #"+ entries[k].index + "  质荷比: " + String.format("%.4f", entries[k].mz) +
                            "   峰强度: " + String.format("%.4f", entries[k].TIC) +
                            "   保留时间: " + String.format("%.4f", entries[k].RT));
                    System.out.println("   #"+ entries[i].index + "  质荷比: " + String.format("%.4f", entries[i].mz) +
                            "   峰强度: " + String.format("%.4f", entries[i].TIC) +
                            "   保留时间: " + String.format("%.4f", entries[i].RT));

                    System.out.println();
                }
            }
        }
    }
}
