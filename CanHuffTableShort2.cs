using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DPCMCompressor {

    public class CanHuffTableShort2 {
        // Return maxCodeLength                            
        public int GetCanHuffmanTable(
            ushort[] inputArray,
            int sizeOfInput,
            ushort[] symbolsOfInterestList,
            int sizeOfSymbolsOfInterest,
            int[] inFrequenceyArray,
            out ushort[] outSymbolArray,
            out int outSizeOfSymbolArray,
            out int[] outCanHuffmanTable,
            out int[] outSymbolVsCodeLength,
            out int[] outCodeLengthVsCodeStart,
            out int[] outSymbolsPerCodeLength,
            bool addDelimiter,
            ushort delimiterSymbol
        ) {
            // Return maxCodeLength
            int[] outFreqArray = null;
            //DataType * outSymbolArray;

            //int numberOfSymbols =
            //    BuildFrequencyTable(
            //        inputArray,
            //        sizeOfInput,
            //        symbolsOfInterestList,
            //        sizeOfSymbolsOfInterest,
            //        out outFreqArray,
            //        out outSymbolArray
            //    );

            int numberOfSymbols = symbolsOfInterestList.Length;
            outFreqArray = inFrequenceyArray;
            outSymbolArray = symbolsOfInterestList;

            // Let us check and add the delimiter
            if (addDelimiter) {
                int[] tempFreqArray = new int[numberOfSymbols + 1];
                ushort[] tempSymbolArray = new ushort[numberOfSymbols + 1];
                Array.Copy(outFreqArray, tempFreqArray, numberOfSymbols);
                //memcpy(tempFreqArray, outFreqArray, (numberOfSymbols) * 4);
                Array.Copy(outSymbolArray, tempSymbolArray, numberOfSymbols);
                //memcpy(tempSymbolArray, outSymbolArray, sizeof(DataType) * (numberOfSymbols));
                tempFreqArray[numberOfSymbols] = 1;
                tempSymbolArray[numberOfSymbols] = delimiterSymbol;
                numberOfSymbols++;
                //delete[] outFreqArray;
                //delete[] outSymbolArray;
                outFreqArray = tempFreqArray;
                outSymbolArray = tempSymbolArray;
            }

            int maxCodeLength =
                CalculateCodeLength(ref outSymbolArray, ref outFreqArray, numberOfSymbols);
            int[] codeLengthPerSymbolArray = outFreqArray;
            int[] symbolsPerCodeLength =
                SymbolsPerCodeLength(codeLengthPerSymbolArray, numberOfSymbols, maxCodeLength);
            int[] symbolStartPerCodeLengthArray = null;
            CalculateSymbolStartForCodeLength(
                codeLengthPerSymbolArray,
                symbolsPerCodeLength,
                maxCodeLength,
                out symbolStartPerCodeLengthArray);

            int[] canHuffmanTable = null;
            ConstructCanHuffmanTable(
                outSymbolArray,
                codeLengthPerSymbolArray,
                numberOfSymbols,
                symbolStartPerCodeLengthArray,
                maxCodeLength,
                out canHuffmanTable);

            //outSizeOfSymbolArray = new int[1];
            outSizeOfSymbolArray = numberOfSymbols;
            outCanHuffmanTable = canHuffmanTable;
            outSymbolVsCodeLength = codeLengthPerSymbolArray;
            outCodeLengthVsCodeStart = symbolStartPerCodeLengthArray;
            outSymbolsPerCodeLength = symbolsPerCodeLength;
            return maxCodeLength;
        }

        private int BuildFrequencyTable(
            ushort[] inputArray,
            int sizeOfInput,
            ushort[] symbolsOfInterestList,
            int sizeOfSymbolsOfInterest,
            out int[] outFreqArray,
            out ushort[] outSymbolArray
        ) {
            bool is16bit = true;// Marshal.SizeOf(typeof(T)) == Marshal.SizeOf(typeof(short));
            int sizeOfArray = is16bit ? (1 << 16) : (1 << 8);
            outSymbolArray = new ushort[sizeOfArray];
            outFreqArray = new int[sizeOfArray];
            int numberOfSymbols = 0;
            for (int i = 0; i < sizeOfInput; i++) {
                ushort currentSymbol = inputArray[i];
                // Check if need to use this symbol
                bool interested = false;
                for (int ij = 0; ij < sizeOfSymbolsOfInterest; ++ij) {
                    // TODO: remove .equals, it is very slow. Instead created specific classes for short and byte
                    if (symbolsOfInterestList[ij] == (currentSymbol)) {
                        interested = true;
                        break;
                    }
                }
                if (interested) {
                    // Check if the current symbol exists in the symbolArray
                    bool symbolExists = false;
                    int existingSymbolPosition = 0;
                    for (int j = 0; j < numberOfSymbols; j++) {
                        if (currentSymbol == (outSymbolArray[j])) {
                            existingSymbolPosition = j;
                            symbolExists = true;
                            break;
                        }
                    }
                    if (symbolExists) {
                        outFreqArray[existingSymbolPosition]++;
                    } else {
                        outSymbolArray[numberOfSymbols] = currentSymbol;
                        outFreqArray[numberOfSymbols] = 1;
                        numberOfSymbols++;
                    }
                }
            }

            // Copy the arrays into exact size arrays
            ushort[] tempSymbolArray = new ushort[numberOfSymbols];
            int[] tempFreqArray = new int[numberOfSymbols];

            // Use memcpy HERE: TODO
            for (int i = 0; i < numberOfSymbols; i++) {
                tempSymbolArray[i] = outSymbolArray[i];
                tempFreqArray[i] = outFreqArray[i];
            }


            outSymbolArray = tempSymbolArray;
            outFreqArray = tempFreqArray;

            return numberOfSymbols;
        }

        void QuickSort(ref ushort[] symbols, ref int[] freq, int leftPos, int rightPos) {
            if (rightPos > leftPos) {
                // Select pivot position.
                // We will select the middle element as the pivot position
                int pivotPos = (rightPos + leftPos) >> 1;

                int newPivotPos = Partition(ref symbols, ref freq, leftPos, rightPos, pivotPos);
                QuickSort(ref symbols, ref freq, leftPos, newPivotPos - 1);
                QuickSort(ref symbols, ref freq, newPivotPos + 1, rightPos);
            }
        }

        int Partition(ref ushort[] symbols, ref int[] freq, int leftPos, int rightPos, int pivotPos) {
            int pivotValue = freq[pivotPos];

            // Swap pivot and right most element
            {
                int rightValue = freq[rightPos];
                freq[rightPos] = pivotValue;
                freq[pivotPos] = rightValue;

                ushort rightSymbol = symbols[rightPos];
                ushort pivotSymbol = symbols[pivotPos];
                symbols[rightPos] = pivotSymbol;
                symbols[pivotPos] = rightSymbol;
            }

            int storePos = leftPos;

            for (int i = leftPos; i < rightPos; ++i) {
                if (freq[i] <= pivotValue) {
                    // Swap 
                    {
                        int iValue = freq[i];
                        int storePosValue = freq[storePos];
                        freq[i] = storePosValue;
                        freq[storePos] = iValue;

                        ushort iSymbol = symbols[i];
                        ushort storePosSymbol = symbols[storePos];
                        symbols[i] = storePosSymbol;
                        symbols[storePos] = iSymbol;
                    }
                    storePos++;
                }
            }

            // Move pivot to its final place
            {
                int rightValue = freq[rightPos];
                int storePosValue = freq[storePos];
                freq[rightPos] = storePosValue;
                freq[storePos] = rightValue;

                ushort rightSymbol = symbols[rightPos];
                ushort storePosSymbol = symbols[storePos];
                symbols[rightPos] = storePosSymbol;
                symbols[storePos] = rightSymbol;
            }
            return storePos;
        }

        int CalculateCodeLength(ref ushort[] symbols, ref int[] freq, int count) {

            // Sort the frequencies in ascending order                              
            QuickSort(ref symbols, ref freq, 0, count - 1);

            // Minimim redudancy code evaluation algorithm written by Alister Moffat and Jyrki Katajainen
            // This code calculates the code lengths in place.
            // http://www.cs.mu.oz.au/~alistair/inplace.c

            int root;                  /* next root node to be used */
            int leaf;                  /* next leaf to be used */
            int next;                  /* next value to be assigned */
            int avbl;                  /* number of available nodes */
            int used;                  /* number of internal nodes */
            int dpth;                  /* current depth of leaves */

            // Check for boundary conditions
            if (count == 0) return 0;
            if (count == 1) {
                // Set the required code lenght as 0   
                freq[0] = 0;
                return 0;
            }

            // First pass
            freq[0] += freq[1]; root = 0; leaf = 2;

            for (next = 1; next < count - 1; ++next) {
                // Select first item for pairing
                if ((leaf >= count) || (freq[root] < freq[leaf])) {
                    freq[next] = freq[root];
                    freq[root++] = next;
                } else {
                    freq[next] = freq[leaf++];
                }
                // Add on the second item
                if ((leaf >= count) || ((root < next) && (freq[root] < freq[leaf]))) {
                    freq[next] += freq[root];
                    freq[root++] = next;
                } else {
                    freq[next] += freq[leaf++];
                }
            }

            // Second pass, right to left
            freq[count - 2] = 0;
            for (next = count - 3; next >= 0; --next) {
                freq[next] = freq[freq[next]] + 1;
            }

            // Third pass, right to left
            avbl = 1;
            used = dpth = 0;
            root = count - 2;
            next = count - 1;
            while (avbl > 0) {
                while ((root >= 0) && (freq[root] == dpth)) {
                    used++;
                    root--;
                }
                while (avbl > used) {
                    freq[next--] = dpth;
                    avbl--;
                }
                avbl = 2 * used;
                dpth++;
                used = 0;
            }
            // Alister Moffat code ends
            return freq[0];
        }

        int[] SymbolsPerCodeLength(int[] codeLengthPerSymbolArray, int numberOfSymbols, int maxCodeLength) {
            int[] numberOfSymbolsArray = new int[maxCodeLength + 1];
            //memset(numberOfSymbolsArray, 0, (maxCodeLength + 1) * 4);
            for (int i = 0; i < numberOfSymbols; ++i) {
                numberOfSymbolsArray[codeLengthPerSymbolArray[i]]++;
            }
            return numberOfSymbolsArray;
        }

        void CalculateSymbolStartForCodeLength(
            int[] codeLengthPerSymbolArray,
            int[] numberOfSymbolsPerCodeLengthArray,
            int maxCodeLength,
            out int[] outSymbolStartPerCodeLengthArray

        ) {
            outSymbolStartPerCodeLengthArray = new int[maxCodeLength + 1];
            //memset(outSymbolStartPerCodeLengthArray, 0, (maxCodeLength + 1) * 4);
            // Wrong comment : The code lengths in 'codeLengthPerSymbolArray' is in decending order
            // It is important that this order is maintained
            int symbolStart = 0;
            int prevCodeLength = 0;
            int numberOfSymbolsForPrevCodeLength = 0;
            for (int i = 1; i < (maxCodeLength + 1); i++) {
                int numberOfSymbols = numberOfSymbolsPerCodeLengthArray[i];
                if (numberOfSymbols != 0) {
                    if (prevCodeLength == 0) {
                        outSymbolStartPerCodeLengthArray[i] = symbolStart;
                        prevCodeLength = i;
                        numberOfSymbolsForPrevCodeLength = numberOfSymbols;
                    } else {
                        outSymbolStartPerCodeLengthArray[i] =
                            (outSymbolStartPerCodeLengthArray[prevCodeLength] +
                            numberOfSymbolsForPrevCodeLength) << (i - prevCodeLength);
                        prevCodeLength = i;
                        numberOfSymbolsForPrevCodeLength = numberOfSymbols;
                    }
                }
            }
        }

        void ConstructCanHuffmanTable(
            ushort[] symbolArray,
            int[] symbolVsCodeLengthArray,
            int numberOfSymbols,
            int[] codeLengthVsCodeStart,
            int maxCodeLength,
            out int[] outCanHuffmanTable

        ) {
            outCanHuffmanTable = new int[numberOfSymbols];
            int[] copyOfCodelengthVsCodeStart = new int[maxCodeLength + 1];
            Array.Copy(codeLengthVsCodeStart, copyOfCodelengthVsCodeStart, maxCodeLength + 1);
            //memcpy(copyOfCodelengthVsCodeStart, codeLengthVsCodeStart, (maxCodeLength + 1) * 4);
            for (int i = 0; i < numberOfSymbols; ++i) {
                outCanHuffmanTable[i] = copyOfCodelengthVsCodeStart[symbolVsCodeLengthArray[i]];
                copyOfCodelengthVsCodeStart[symbolVsCodeLengthArray[i]]++;
            }
        }
    }
}
