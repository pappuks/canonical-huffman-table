using DPCMCompressor;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DeltaRLEHuffCompressor {
    public class HuffmanCompressor {

        const int HUFFMAN_SYMBOLS = 500;
        const int DIV_FACTOR = 200;

        public static byte[] Compress(ushort[] input, int pixelDepth) {
            ushort sections = input[0];

            int headerSize = 1 /* sections */+ 2 /* total length 32 bit*/ + 2 /* section length 32 bit */ +
                sections * 2 /*delta section sizes*/ + sections * 2 /*RLE section sizes*/ + sections * 2 /* Huffman section sizes*/;

            // convert header size to bytes
            int headerSizeBytes = headerSize * 2;

            int Regions = (1 << pixelDepth);

            int delimiterForCompressDecompress = (int)((1 << (pixelDepth)) - 1);

            ushort[] shortSymbolsOfInterest = null;
            int[] freqOfSymbols = null;

            GenerateFrequencies(input, headerSize, Regions, delimiterForCompressDecompress, out shortSymbolsOfInterest, out freqOfSymbols);

            int outSizeOfSymbolArray = 0;
            int[] outCanHuffmanTable = null;
            int[] outSymbolVsCodeLength = null;
            int[] outCodeLengthVsCodeStart = null;
            int[] outSymbolsPerCodeLength = null;
            int maxCodeLength = 0;

            int indexOfDelimiter = 0;
            ushort[] outSymbolArrayShort = null;
            CanHuffTableShort2 canHuffmanTable = new CanHuffTableShort2();
            maxCodeLength =
                canHuffmanTable.GetCanHuffmanTable(
                    input,
                    input.Length,
                    shortSymbolsOfInterest,
                    shortSymbolsOfInterest.Length,
                    freqOfSymbols,
                    out outSymbolArrayShort,
                    out outSizeOfSymbolArray,
                    out outCanHuffmanTable,
                    out outSymbolVsCodeLength,
                    out outCodeLengthVsCodeStart,
                    out outSymbolsPerCodeLength,
                    false, // Don't add delimiter as it is already part of the symbol list
                    (ushort)delimiterForCompressDecompress);

            for (int i = 0; i < outSizeOfSymbolArray; ++i) {
                if (
                    (ushort)delimiterForCompressDecompress ==
                    (ushort)outSymbolArrayShort[i]
                ) {
                    indexOfDelimiter = i;
                    break;
                }
            }

            byte[] outputDataPtr = new byte[input.Length * 4 + headerSizeBytes];

            // Copy the header
            for(int i = 0; i < headerSize; i++) {
                outputDataPtr[2 * i] = (byte)input[i];
                outputDataPtr[(2 * i) + 1] = (byte)(input[i] >> 8);
            }

            int outDataCounter = headerSizeBytes;

            int imageSize = input.Length - headerSize;

            // Write details to output array

            // Store the uncompressed image size
            uint imageSizeCopy = (uint)imageSize * 2/*imageSize*/;
            outputDataPtr[outDataCounter + 3] = (byte)(imageSizeCopy >> 24);
            outputDataPtr[outDataCounter + 2] = (byte)(imageSizeCopy >> 16);
            outputDataPtr[outDataCounter + 1] = (byte)(imageSizeCopy >> 8);
            outputDataPtr[outDataCounter + 0] = (byte)(imageSizeCopy);
            //            memcpy(outputDataPtr + outDataCounter,&imageSizeCopy, 4);
            outDataCounter += 4;
            // Store the delta array size
            outputDataPtr[outDataCounter + 3] = (byte)((uint)imageSize >> 24);
            outputDataPtr[outDataCounter + 2] = (byte)((uint)imageSize >> 16);
            outputDataPtr[outDataCounter + 1] = (byte)((uint)imageSize >> 8);
            outputDataPtr[outDataCounter + 0] = (byte)((uint)imageSize);
            //memcpy(outputDataPtr + outDataCounter,&deltaArrayCounter, 4);
            outDataCounter += 4;
            uint pixelDepthCopy = (uint)pixelDepth;
            // Store the pixelDepth
            outputDataPtr[outDataCounter + 3] = (byte)((uint)pixelDepthCopy >> 24);
            outputDataPtr[outDataCounter + 2] = (byte)((uint)pixelDepthCopy >> 16);
            outputDataPtr[outDataCounter + 1] = (byte)((uint)pixelDepthCopy >> 8);
            outputDataPtr[outDataCounter + 0] = (byte)((uint)pixelDepthCopy);
            //memcpy(outputDataPtr + outDataCounter,&pixelDepthCopy, 4);
            outDataCounter += 4;
            // Store the size of symbols array
            outputDataPtr[outDataCounter + 3] = (byte)((uint)outSizeOfSymbolArray >> 24);
            outputDataPtr[outDataCounter + 2] = (byte)((uint)outSizeOfSymbolArray >> 16);
            outputDataPtr[outDataCounter + 1] = (byte)((uint)outSizeOfSymbolArray >> 8);
            outputDataPtr[outDataCounter + 0] = (byte)((uint)outSizeOfSymbolArray);
            outDataCounter += 4;
            // Store the symbol list            
            for (int i = 0; i < outSizeOfSymbolArray; i++) {
                outputDataPtr[outDataCounter + (2 * i) + 1] = (byte)(outSymbolArrayShort[i] >> 8);
                outputDataPtr[outDataCounter + (2 * i)] = (byte)(outSymbolArrayShort[i]);
            }
            outDataCounter += outSizeOfSymbolArray * 2;
            // Store the huffman table
            for (int i = 0; i < outSizeOfSymbolArray; i++) {
                outputDataPtr[outDataCounter + 3 + (4 * i)] = (byte)((uint)outCanHuffmanTable[i] >> 24);
                outputDataPtr[outDataCounter + 2 + (4 * i)] = (byte)((uint)outCanHuffmanTable[i] >> 16);
                outputDataPtr[outDataCounter + 1 + (4 * i)] = (byte)((uint)outCanHuffmanTable[i] >> 8);
                outputDataPtr[outDataCounter + 0 + (4 * i)] = (byte)((uint)outCanHuffmanTable[i]);
            }
            outDataCounter += outSizeOfSymbolArray * 4;
            // Store symbol vs code length
            for (int i = 0; i < outSizeOfSymbolArray; i++) {
                outputDataPtr[outDataCounter + 3 + (4 * i)] = (byte)((uint)outSymbolVsCodeLength[i] >> 24);
                outputDataPtr[outDataCounter + 2 + (4 * i)] = (byte)((uint)outSymbolVsCodeLength[i] >> 16);
                outputDataPtr[outDataCounter + 1 + (4 * i)] = (byte)((uint)outSymbolVsCodeLength[i] >> 8);
                outputDataPtr[outDataCounter + 0 + (4 * i)] = (byte)((uint)outSymbolVsCodeLength[i]);
            }
            outDataCounter += outSizeOfSymbolArray * 4;
            // Store max code length
            outputDataPtr[outDataCounter + 3] = (byte)((uint)maxCodeLength >> 24);
            outputDataPtr[outDataCounter + 2] = (byte)((uint)maxCodeLength >> 16);
            outputDataPtr[outDataCounter + 1] = (byte)((uint)maxCodeLength >> 8);
            outputDataPtr[outDataCounter + 0] = (byte)((uint)maxCodeLength);
            outDataCounter += 4;
            // Store code length vs code start
            // The size of the code length to code start array is maxCodeLength + 1
            for (int i = 0; i < (maxCodeLength + 1); i++) {
                outputDataPtr[outDataCounter + 3 + (4 * i)] = (byte)((uint)outCodeLengthVsCodeStart[i] >> 24);
                outputDataPtr[outDataCounter + 2 + (4 * i)] = (byte)((uint)outCodeLengthVsCodeStart[i] >> 16);
                outputDataPtr[outDataCounter + 1 + (4 * i)] = (byte)((uint)outCodeLengthVsCodeStart[i] >> 8);
                outputDataPtr[outDataCounter + 0 + (4 * i)] = (byte)((uint)outCodeLengthVsCodeStart[i]);
            }
            outDataCounter += (maxCodeLength + 1) * 4;
            // Store symbols per code length
            for (int i = 0; i < (maxCodeLength + 1); i++) {
                outputDataPtr[outDataCounter + 3 + (4 * i)] = (byte)((uint)outSymbolsPerCodeLength[i] >> 24);
                outputDataPtr[outDataCounter + 2 + (4 * i)] = (byte)((uint)outSymbolsPerCodeLength[i] >> 16);
                outputDataPtr[outDataCounter + 1 + (4 * i)] = (byte)((uint)outSymbolsPerCodeLength[i] >> 8);
                outputDataPtr[outDataCounter + 0 + (4 * i)] = (byte)((uint)outSymbolsPerCodeLength[i]);
            }
            outDataCounter += (maxCodeLength + 1) * 4;

            int readSectionLength = headerSize;
            
            for (int sec = 0; sec < sections; sec++) {
                int sectionEnding = input[(2 * sec) + 5 + (2 * sections)] + (input[(2 * sec) + 1 + 5 + (2 * sections)] << 16);
                int startIndex = readSectionLength;
                int endIndex = sectionEnding;
                readSectionLength = sectionEnding;

                int bitsOccupied = 0;
                byte charToWrite = 0;
                for (int i = startIndex; i < endIndex; ++i) {
                    ushort currentSymbol = 0;
                    currentSymbol = input[i];

                    bool appendDelimiter = true;
                    int indexOfSymbol = 0; // TODO: Here default this to indexOfDelimiter

                    if (currentSymbol == (ushort)delimiterForCompressDecompress) {
                        appendDelimiter = true;
                        indexOfSymbol = indexOfDelimiter;
                    } else {
                        for (int j = (outSizeOfSymbolArray - 1); j >= 0; j--) {
                            if (currentSymbol == (ushort)(outSymbolArrayShort[j])) {
                                appendDelimiter = false;
                                indexOfSymbol = j;
                                break;
                            }
                        }
                    }

                    int codeLength = appendDelimiter ? outSymbolVsCodeLength[indexOfDelimiter] : outSymbolVsCodeLength[indexOfSymbol];
                    int code = appendDelimiter ? outCanHuffmanTable[indexOfDelimiter] : outCanHuffmanTable[indexOfSymbol];

                    // Write the code to buffer;
                    int bitsLeft = 8 - bitsOccupied;
                    if (bitsLeft > codeLength) {
                        charToWrite =
                            (byte)(charToWrite | (byte)(((byte)code) << (8 - bitsOccupied - codeLength)));
                        bitsOccupied += codeLength;
                    } else {
                        int leftCodeLength = codeLength;
                        while (leftCodeLength > 0) {
                            if (bitsLeft <= leftCodeLength) {
                                charToWrite =
                                    (byte)(charToWrite | (byte)(code >> (leftCodeLength - bitsLeft)));
                                leftCodeLength -= bitsLeft;
                                bitsOccupied += bitsLeft;
                                if (bitsOccupied == 8) {
                                    outputDataPtr[outDataCounter++] = charToWrite;
                                    charToWrite = 0;
                                    bitsOccupied = 0;
                                    bitsLeft = 8;
                                }
                            } else {
                                charToWrite =
                                    (byte)(charToWrite | (byte)(((byte)code) << (8 - bitsOccupied - leftCodeLength)));
                                bitsOccupied += leftCodeLength;
                                leftCodeLength = 0; // To break the loop
                            }
                        }
                    }

                    if (appendDelimiter) {
                        // Now write the symbol
                        bitsLeft = 8 - bitsOccupied;
                        int symbolLength = pixelDepth; // WE HAVE TO DOUBLE CHECK THIS ASSIGNMENT(pixelSize > 8) ? 16 : 8;
                        if (bitsLeft > symbolLength) {
                            charToWrite =
                                (byte)(charToWrite | (byte)(((byte)currentSymbol) << (8 - bitsOccupied - symbolLength)));
                            bitsOccupied += symbolLength;
                        } else {
                            int leftCodeLength = symbolLength;
                            while (leftCodeLength > 0) {
                                if (bitsLeft <= leftCodeLength) {
                                    charToWrite =
                                        (byte)(charToWrite | (byte)(currentSymbol >> (leftCodeLength - bitsLeft)));
                                    leftCodeLength -= bitsLeft;
                                    bitsOccupied += bitsLeft;
                                    if (bitsOccupied == 8) {
                                        outputDataPtr[outDataCounter++] = charToWrite;
                                        charToWrite = 0;
                                        bitsOccupied = 0;
                                        bitsLeft = 8;
                                    }
                                } else {
                                    charToWrite =
                                        (byte)(charToWrite | (byte)(((byte)currentSymbol) << (8 - bitsOccupied - leftCodeLength)));
                                    bitsOccupied += leftCodeLength;
                                    leftCodeLength = 0; // To break the loop
                                }
                            }
                        }
                    }
                }

                // Flush the last byte
                outputDataPtr[outDataCounter++] = charToWrite;

                // Store the size
                outputDataPtr[(4 * sec) + 10 + (8 * sections)] = (byte)(outDataCounter);
                outputDataPtr[(4 * sec) + 1 + 10 + (8 * sections)] = (byte)((outDataCounter) >> 8);
                outputDataPtr[(4 * sec) + 2 + 10 + (8 * sections)] = (byte)((outDataCounter) >> 16);
                outputDataPtr[(4 * sec) + 3 + 10 + (8 * sections)] = (byte)((outDataCounter) >> 24);
            }

            byte[] retArray = new byte[outDataCounter];
            Array.Copy(outputDataPtr, retArray, outDataCounter);
            return retArray;
        }

        private static void GenerateFrequencies(
            ushort[] input, 
            int headerSize, 
            int Regions, 
            int delimiterForCompressDecompress,
            out ushort[] shortSymbolsOfInterest,
            out int[] freqOfSymbols
        ) {
            int[] distributionArray = new int[Regions];

            int maxRegionValue = 0;
            int maxRegion = 0;
            for (int i = headerSize; i < input.Length; i++) {
                distributionArray[input[i]]++;
                if (distributionArray[input[i]] > maxRegionValue) {
                    maxRegionValue = distributionArray[input[i]];
                    maxRegion = input[i];
                }
            }

            List<KeyValuePair<ushort, int>> shortSymbolsOfInterestList = new List<KeyValuePair<ushort, int>>();

            for (int i = 0; i < Regions; i++) {
                if (distributionArray[i] > 0) {
                    if (i != delimiterForCompressDecompress) {
                        shortSymbolsOfInterestList.Add(new KeyValuePair<ushort, int>((ushort)i, distributionArray[i]));
                    }
                }
            }

            shortSymbolsOfInterestList.Sort((x, y) => { return y.Value - x.Value; }); // Sort in descending order

            // Take only symbols which fall within 1/100 of the max freq symbol.
            int maxFrequency = shortSymbolsOfInterestList[0].Value;

            int minAllowedFrequency = maxFrequency / DIV_FACTOR;

            for (int i = 0; i < shortSymbolsOfInterestList.Count; i++) {
                if (shortSymbolsOfInterestList[i].Value < minAllowedFrequency) {
                    shortSymbolsOfInterestList.RemoveRange(i, shortSymbolsOfInterestList.Count - i);
                    break;
                }
            }

            // Take only the first 500 symbols
            if (shortSymbolsOfInterestList.Count > HUFFMAN_SYMBOLS) {
                shortSymbolsOfInterestList.RemoveRange(HUFFMAN_SYMBOLS, shortSymbolsOfInterestList.Count - HUFFMAN_SYMBOLS);
            }

            // Add the delimiter correctly -- START
            long selectedSymbolCount = 0;
            for (int i = 0; i < shortSymbolsOfInterestList.Count; i++) {
                selectedSymbolCount += shortSymbolsOfInterestList[i].Value;
            }

            int delimiterCount = (int)((input.Length - headerSize) - selectedSymbolCount);
            shortSymbolsOfInterestList.Add(new KeyValuePair<ushort, int>((ushort)delimiterForCompressDecompress, (int)delimiterCount));

            // Sort again after adding delimiter symbol
            shortSymbolsOfInterestList.Sort((x, y) => { return y.Value - x.Value; }); // Sort in descending order

            // Add the delimiter correctly -- END

            shortSymbolsOfInterest = new ushort[shortSymbolsOfInterestList.Count];
            freqOfSymbols = new int[shortSymbolsOfInterestList.Count];

            for (int i = 0; i < shortSymbolsOfInterestList.Count; i++) {
                shortSymbolsOfInterest[i] = shortSymbolsOfInterestList[i].Key;
                freqOfSymbols[i] = shortSymbolsOfInterestList[i].Value;
            }            
            
        }

        private int sections;
        private int headerSizeInBytes;
        private int pixelDepth;
        private int sizeOfSymbolsArray;
        ushort[] shortSymbolArray;
        int[] canHuffmanTable;
        int[] symbolVsCodeLength;
        private int maxCodeLength;
        int[] codeLengthVsCodeStart;
        int[] symbolsPerCodeLength;
        int[] leftShiftedCodeStartArray;
        int[][] codeLengthToSymbolTable;
        int[] correctCodeLengthsArray;
        int[] maskForByteCodeLengths;
        int delimiterForCompressDecompress;
        ushort deltaThreshold;
        ushort delimiterForCodecShort;

        public int ReadDecodeInfo(byte[] input) {
            const int sizeOfMaxCodeLengthArray = 32;

            sections = input[0] + (input[1] << 8);

            int headerSize = 1 /* sections */+ 2 /* total length 32 bit*/ + 2 /* section length 32 bit */ +
                sections * 2 /*delta section sizes*/ + sections * 2 /*RLE section sizes*/ + sections * 2 /* Huffman section sizes*/;

            // convert header size to bytes
            headerSizeInBytes = headerSize * 2;


            // Retrieve the stored information from the image
            byte[] inputDataPtrChar = input;
            int inputDataCounter = headerSizeInBytes;

            // Read image size
            int imageSize =
                (int)(((uint)inputDataPtrChar[inputDataCounter + 3] << 24) +
                ((uint)inputDataPtrChar[inputDataCounter + 2] << 16) +
                ((uint)inputDataPtrChar[inputDataCounter + 1] << 8) +
                inputDataPtrChar[inputDataCounter + 0]);

            inputDataCounter += 4;
            // Read delta array size
            int deltaArraySize =
                (int)(((uint)inputDataPtrChar[inputDataCounter + 3] << 24) +
                ((uint)inputDataPtrChar[inputDataCounter + 2] << 16) +
                ((uint)inputDataPtrChar[inputDataCounter + 1] << 8) +
                inputDataPtrChar[inputDataCounter + 0]);
            inputDataCounter += 4;
            // Read pixelDepth
            pixelDepth =
                (int)(((uint)inputDataPtrChar[inputDataCounter + 3] << 24) +
                ((uint)inputDataPtrChar[inputDataCounter + 2] << 16) +
                ((uint)inputDataPtrChar[inputDataCounter + 1] << 8) +
                inputDataPtrChar[inputDataCounter + 0]);

            inputDataCounter += 4;

            // Read the size of symbols array
            sizeOfSymbolsArray =
                (int)(((uint)inputDataPtrChar[inputDataCounter + 3] << 24) +
                ((uint)inputDataPtrChar[inputDataCounter + 2] << 16) +
                ((uint)inputDataPtrChar[inputDataCounter + 1] << 8) +
                inputDataPtrChar[inputDataCounter + 0]);

            inputDataCounter += 4;
            // Read the symbols array
            int sizeOfActualSymbolArray = (sizeOfSymbolsArray * 2);

            shortSymbolArray = new ushort[sizeOfSymbolsArray];

            for (int i = 0; i < sizeOfSymbolsArray; i++) {
                shortSymbolArray[i] = (ushort)
                    (((ushort)inputDataPtrChar[inputDataCounter + (2 * i) + 1] << 8) +
                    (ushort)inputDataPtrChar[inputDataCounter + (2 * i)]);
            }
            inputDataCounter += sizeOfActualSymbolArray;
            // Read the huffman table
            canHuffmanTable = new int[sizeOfSymbolsArray];
            for (int i = 0; i < sizeOfSymbolsArray; i++) {
                canHuffmanTable[i] =
                (int)(((uint)inputDataPtrChar[inputDataCounter + 3 + (4 * i)] << 24) +
                ((uint)inputDataPtrChar[inputDataCounter + 2 + (4 * i)] << 16) +
                ((uint)inputDataPtrChar[inputDataCounter + 1 + (4 * i)] << 8) +
                inputDataPtrChar[inputDataCounter + 0 + (4 * i)]);
            }
            inputDataCounter += sizeOfSymbolsArray * 4;
            // Read symbol vs code length
            symbolVsCodeLength = new int[sizeOfSymbolsArray];
            for (int i = 0; i < sizeOfSymbolsArray; i++) {
                symbolVsCodeLength[i] =
                (int)(((uint)inputDataPtrChar[inputDataCounter + 3 + (4 * i)] << 24) +
                ((uint)inputDataPtrChar[inputDataCounter + 2 + (4 * i)] << 16) +
                ((uint)inputDataPtrChar[inputDataCounter + 1 + (4 * i)] << 8) +
                inputDataPtrChar[inputDataCounter + 0 + (4 * i)]);
            }
            inputDataCounter += sizeOfSymbolsArray * 4;
            // Read maxCode Length
            maxCodeLength =
                (int)(((uint)inputDataPtrChar[inputDataCounter + 3] << 24) +
                ((uint)inputDataPtrChar[inputDataCounter + 2] << 16) +
                ((uint)inputDataPtrChar[inputDataCounter + 1] << 8) +
                inputDataPtrChar[inputDataCounter + 0]);
            inputDataCounter += 4;
            // Read code length vs code start
            codeLengthVsCodeStart = new int[sizeOfMaxCodeLengthArray];
            for (int i = 0; i < (maxCodeLength + 1); i++) {
                codeLengthVsCodeStart[i] =
                (int)(((uint)inputDataPtrChar[inputDataCounter + 3 + (4 * i)] << 24) +
                ((uint)inputDataPtrChar[inputDataCounter + 2 + (4 * i)] << 16) +
                ((uint)inputDataPtrChar[inputDataCounter + 1 + (4 * i)] << 8) +
                inputDataPtrChar[inputDataCounter + 0 + (4 * i)]);
            }
            inputDataCounter += (maxCodeLength + 1) * 4;
            //Read symbols per code length
            symbolsPerCodeLength = new int[sizeOfMaxCodeLengthArray];
            for (int i = 0; i < (maxCodeLength + 1); i++) {
                symbolsPerCodeLength[i] =
                (int)(((uint)inputDataPtrChar[inputDataCounter + 3 + (4 * i)] << 24) +
                ((uint)inputDataPtrChar[inputDataCounter + 2 + (4 * i)] << 16) +
                ((uint)inputDataPtrChar[inputDataCounter + 1 + (4 * i)] << 8) +
                inputDataPtrChar[inputDataCounter + 0 + (4 * i)]);
            }
            inputDataCounter += (maxCodeLength + 1) * 4;

            leftShiftedCodeStartArray = new int[sizeOfMaxCodeLengthArray];

            for (int i = 1; i < (maxCodeLength + 1); ++i) {
                leftShiftedCodeStartArray[i] = codeLengthVsCodeStart[i] << (maxCodeLength - i);
            }
            // Set one final value greater than the max left shifted code start value to act as delimiter
            leftShiftedCodeStartArray[maxCodeLength + 1] = leftShiftedCodeStartArray[maxCodeLength] + 1 + symbolsPerCodeLength[maxCodeLength];

            //Prepare a code length table which will point to arrays which contain the symbols
            codeLengthToSymbolTable = new int[maxCodeLength + 1][];

            // From the canHuffmanTable  and code length table we have to prepare a per codelength symbol table
            correctCodeLengthsArray = new int[sizeOfMaxCodeLengthArray];
            int tempCounter = 0;
            // Initialize the codeLength to symbol table
            for (int i = 1; i < (maxCodeLength + 1); ++i) {
                int symbolsForThisCodeLength = symbolsPerCodeLength[i];
                int[] tempSymbolList = new int[symbolsForThisCodeLength];
                codeLengthToSymbolTable[i] = tempSymbolList;
                if (symbolsForThisCodeLength > 0) {
                    // Qualified codeLength
                    correctCodeLengthsArray[tempCounter++] = i;
                }
            }
            // Fill a large maxCodeLength + 1 value for delimiter purpose
            correctCodeLengthsArray[tempCounter++] = maxCodeLength + 1;

            int minCodeLength = maxCodeLength;

            // Populate the codeLength to symbol table    
            for (int j = 0; j < sizeOfSymbolsArray; ++j) {
                int currentCodeLength = symbolVsCodeLength[j];
                minCodeLength = (currentCodeLength > minCodeLength) ? minCodeLength : currentCodeLength;
                int[] tempSymbolList = codeLengthToSymbolTable[currentCodeLength];
                ushort foundSymbol = 0;
                foundSymbol = (ushort)shortSymbolArray[j];
                int codeForSymbol = canHuffmanTable[j];
                int codeStartForLength = codeLengthVsCodeStart[currentCodeLength];
                tempSymbolList[codeForSymbol - codeStartForLength] = (ushort)foundSymbol;
            }

            maskForByteCodeLengths = new int[33];
            for (int k = 0; k < 33; ++k) {
                maskForByteCodeLengths[k] = (int)(0xffffffff >> (32 - k));
            }

            // Start decoding
            // At a time we read 'maxCodeLength' bits from the buffer and compare that into the leftShiftedCodeStart array.
            // If the value is less than the higher row then we use that index to find the code length and the use that to find the symbol.

            delimiterForCompressDecompress = (int)((1 << (pixelDepth)) - 1);
            deltaThreshold = (ushort)((1 << (pixelDepth - 1)) - 1);
            delimiterForCodecShort = (ushort)delimiterForCompressDecompress;            

            ushort[] shortDecompDeltaBuffer = new ushort[deltaArraySize];

            ushort[] decompDeltaBufferShort = shortDecompDeltaBuffer;

            return inputDataCounter;
        }        

        public void Decompress(byte[] inputDataPtrChar, int section, int firstSectionStart, ushort[] decompDeltaBufferShort) {
            
            int readStart = firstSectionStart;

            // Read the Huffman section start and end
            if (section > 0) {
                readStart = inputDataPtrChar[(4 * (section - 1)) + 10 + (8 * sections)] +
                    (inputDataPtrChar[(4 * (section - 1)) + 1 + 10 + (8 * sections)] << 8) +
                    (inputDataPtrChar[(4 * (section - 1)) + 2 + 10 + (8 * sections)] << 16) +
                    (inputDataPtrChar[(4 * (section - 1)) + 3 + 10 + (8 * sections)] << 24);
            }
            int inputDataCounter = readStart;

            int readEnd = inputDataPtrChar[(4 * section) + 10 + (8 * sections)] +
                (inputDataPtrChar[(4 * section) + 1 + 10 + (8 * sections)] << 8) +
                (inputDataPtrChar[(4 * section) + 2 + 10 + (8 * sections)] << 16) +
                (inputDataPtrChar[(4 * section) + 3 + 10 + (8 * sections)] << 24);

            int writeStart = headerSizeInBytes / 2;

            // Read RLE section start and end indices
            if (section > 0) {
                writeStart = inputDataPtrChar[(4 * (section - 1)) + 10 + (4 * sections)] +
                    (inputDataPtrChar[(4 * (section - 1)) + 1 + 10 + (4 * sections)] << 8) +
                    (inputDataPtrChar[(4 * (section - 1)) + 2 + 10 + (4 * sections)] << 16) +
                    (inputDataPtrChar[(4 * (section - 1)) + 3 + 10 + (4 * sections)] << 24);
            }
            int writeEnd = inputDataPtrChar[(4 * section) + 10 + (4 * sections)] +
                (inputDataPtrChar[(4 * section) + 1 + 10 + (4 * sections)] << 8) +
                (inputDataPtrChar[(4 * section) + 2 + 10 + (4 * sections)] << 16) +
                (inputDataPtrChar[(4 * section) + 3 + 10 + (4 * sections)] << 24);

            int decompBufferCounter = writeStart;
            uint decodeReadInt = 0;
            int decodeByteBitsLeft = 0;
            int readMaxCodeLengthBits = 0;
            int symbolRawLength = pixelDepth;
            int maxCodeLengthMask = (int)(0xffffffff >> (32 - maxCodeLength));
            int bitsLeftToFillInMaxCodeLengthBuffer = maxCodeLength;
            bool wasPreviousSymbolDelimiter = false;
            //
            int doubleDecodeIntBitsLeft = 0;
            uint doubleDecodeReadInt = 0;
            //
            while (decompBufferCounter < writeEnd) {                
                if (decodeByteBitsLeft < bitsLeftToFillInMaxCodeLengthBuffer) {
                    // Left allign the available bits
                    if (decodeByteBitsLeft == 0) {
                        decodeReadInt = 0;
                    }
                    decodeReadInt = (decodeReadInt << (32 - decodeByteBitsLeft));
                    while (decodeByteBitsLeft != 32) {
                        if (doubleDecodeIntBitsLeft == 0) {
                            if ((inputDataCounter + 3) < readEnd) {
                                doubleDecodeReadInt =
                                    (uint)((uint)inputDataPtrChar[inputDataCounter] << 24) +
                                    ((uint)inputDataPtrChar[inputDataCounter + 1] << 16) +
                                    ((uint)inputDataPtrChar[inputDataCounter + 2] << 8) +
                                    (inputDataPtrChar[inputDataCounter + 3]);
                                inputDataCounter += 4;
                            } else {
                                if ((inputDataCounter + 2) < readEnd) {
                                    doubleDecodeReadInt =
                                        (uint)((uint)inputDataPtrChar[inputDataCounter] << 24) +
                                        ((uint)inputDataPtrChar[inputDataCounter + 1] << 16) +
                                        ((uint)inputDataPtrChar[inputDataCounter + 2] << 8);
                                    inputDataCounter += 3;
                                } else if ((inputDataCounter + 1) < readEnd) {
                                    doubleDecodeReadInt =
                                        (uint)((uint)inputDataPtrChar[inputDataCounter] << 24) +
                                        ((uint)inputDataPtrChar[inputDataCounter + 1] << 16);
                                    inputDataCounter += 2;
                                } else if ((inputDataCounter) < readEnd) {
                                    doubleDecodeReadInt =
                                        (uint)((uint)inputDataPtrChar[inputDataCounter] << 24);
                                    inputDataCounter += 1;
                                } else {
                                    // Nothing to read from the input buffer, just put 0 and let the bits left be
                                    // equal to 32.
                                    doubleDecodeReadInt = 0;
                                }
                            }

                            doubleDecodeIntBitsLeft = 32;
                        }
                        int temp32MinusVal = (32 - decodeByteBitsLeft);

                        int leftShiftVal = temp32MinusVal - doubleDecodeIntBitsLeft;
                        if (leftShiftVal >= 0) {
                            decodeReadInt = decodeReadInt | (uint)((uint)(doubleDecodeReadInt << leftShiftVal) & maskForByteCodeLengths[temp32MinusVal]);
                            decodeByteBitsLeft += doubleDecodeIntBitsLeft;
                            doubleDecodeIntBitsLeft = 0;
                        } else {
                            int rightShiftVal = -leftShiftVal;
                            decodeReadInt = decodeReadInt | (uint)((doubleDecodeReadInt >> rightShiftVal) & maskForByteCodeLengths[temp32MinusVal]);
                            doubleDecodeIntBitsLeft -= temp32MinusVal;
                            decodeByteBitsLeft = 32;
                            break;
                        }
                    }
                }                
                readMaxCodeLengthBits = (int)(decodeReadInt >> (decodeByteBitsLeft - bitsLeftToFillInMaxCodeLengthBuffer)) & maskForByteCodeLengths[bitsLeftToFillInMaxCodeLengthBuffer];

                if (!wasPreviousSymbolDelimiter) {
                    int tempCodeLengthCounter = 0;
                    int thisCodeLength = correctCodeLengthsArray[tempCodeLengthCounter];
                    while (true) {
                        tempCodeLengthCounter++;
                        int nextCodeLength = correctCodeLengthsArray[tempCodeLengthCounter];
                        if (
                            (leftShiftedCodeStartArray[nextCodeLength] > readMaxCodeLengthBits)
                        ) {
                            // Set this to minCodeLength as default. We will fill the correct code length below.
                            int correctCodeLength = thisCodeLength;

                            int foundSymbolCode = readMaxCodeLengthBits >> (maxCodeLength - correctCodeLength);
                            int[] tempSymbolTable = codeLengthToSymbolTable[correctCodeLength];
                            int tempCodeStart = codeLengthVsCodeStart[correctCodeLength];
                            ushort foundSymbol = (ushort)tempSymbolTable[foundSymbolCode - tempCodeStart];
                            
                            if (foundSymbol != delimiterForCodecShort) {
                                decompDeltaBufferShort[decompBufferCounter++] = foundSymbol; // Already ushort
                                bitsLeftToFillInMaxCodeLengthBuffer = maxCodeLength;
                                decodeByteBitsLeft -= correctCodeLength;
                                wasPreviousSymbolDelimiter = false;
                            } else {
                                if ((maxCodeLength - correctCodeLength) >= symbolRawLength) {
                                    // We can read the symbol from the existing bits in readMaxCodeLengthBits
                                    foundSymbol = (ushort)((readMaxCodeLengthBits >> (maxCodeLength - correctCodeLength - symbolRawLength)) & (maskForByteCodeLengths[symbolRawLength]));
                                    decompDeltaBufferShort[decompBufferCounter++] = foundSymbol; // Already ushort
                                    bitsLeftToFillInMaxCodeLengthBuffer = maxCodeLength;
                                    decodeByteBitsLeft -= (correctCodeLength + symbolRawLength);
                                    wasPreviousSymbolDelimiter = false;
                                } else {
                                    // We need to read the symbol from the stream
                                    wasPreviousSymbolDelimiter = true;                                    
                                    bitsLeftToFillInMaxCodeLengthBuffer = symbolRawLength;
                                    decodeByteBitsLeft -= correctCodeLength;                                    
                                }
                            }
                            break;
                        } else {
                            thisCodeLength = nextCodeLength;
                        }
                    }
                } else {
                    // Our last symbol was a delimiter, so the value of readMaxCodeLengthBits is the actual symbol
                    decompDeltaBufferShort[decompBufferCounter++] = (ushort)readMaxCodeLengthBits;                    
                    bitsLeftToFillInMaxCodeLengthBuffer = maxCodeLength;
                    decodeByteBitsLeft -= symbolRawLength;
                    readMaxCodeLengthBits = 0;
                    wasPreviousSymbolDelimiter = false;
                }
            }            
        }
    }
}
