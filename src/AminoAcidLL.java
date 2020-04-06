class AminoAcidLL{
  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;

  AminoAcidLL(){

  }


  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon 
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon){
    //method call to get the amino acid character of an specific codon
    this.aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
    //method call to get the possible codons [array] of the amino acid
    this.codons = AminoAcidResources.getCodonListForAminoAcid(this.aminoAcid);
    //initializes the counts for each codon
    this.counts = new int[codons.length];

    Counter(inCodon);

    this.next = null;
  
  }

  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops, 
   * if not passes the task to the next node. 
   * If there is no next node, add a new node to the list that would contain the codon. 
   */
  private void addCodon(String inCodon){

    //This will compare this amino acid with the result from the codons already passed
    //and if its the same it will add to the counter
    if (aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)) {
      this.Counter(inCodon);
    }
    //If the first check did not pass this will check the other nodes (if there are any)
    else if (next != null) {
      next.addCodon(inCodon);
    }
    // If we never found a match we will now create a new node
    else {
      this.next = new AminoAcidLL(inCodon);
    }
  }
  /**************************************HELPER METHOD******************************************/
  //This will count the codon usage
  private void Counter(String inCodon) {
    for (int i = 0; i < this.codons.length; i++) {
      if (this.codons[i].equals(inCodon)) {
        this.counts[i]++;
      }
    }
  }


  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
    int sum = 0;
    for (int i = 0; i < counts.length; i++) {
      sum += counts[i];
    }
    return sum;

}

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    int difference = 0;
    for (int i = 0; i < codons.length; i++) {
      difference += Math.abs(counts[i] - inList.counts[i]);
    }
    return difference;

  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts. 
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){
    if(!inList.isSorted())
      return -1;
    int difference = 88;
    if(inList == null) {
      difference += totalCount();
      if(next != null)
        difference += next.aminoAcidCompare(inList);
    }
    if(next == null) {
      difference += inList.totalCount();
      if(inList.next != null)
        difference += aminoAcidCompare(inList.next);
    }
    difference = totalCount();
    return difference;
  }

  /********************************************************************************************/
  /* Same ad above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList){
    AminoAcidLL sortedInList = sort(inList); //sorts inList as method description indicates
    int difference = codonDiff(inList);
    if (next == null) {
      return difference;
    } else {
      difference = difference + next.codonCompare(inList.next);
    }
    return difference;
  }


  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){
    //base case
    if (next == null) {
      return new char[]{aminoAcid};
    }

    //recursive call to "append" the aminoacids
    char[] temp = next.aminoAcidList();
    //array that stores each of the amino acid characters
    char[] give = new char[temp.length+1];

    //copies temp values in ret
    give[0] = aminoAcid;
    for(int i = 0; i < temp.length; i++){
      give[i+1] = temp[i];
    }

    return give;
  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){
    //base case
    if(next == null){
      return new int[]{this.totalCount()};
    }

    //recursive call to get the total counts of the aminoacid
    int[] temp = next.aminoAcidCounts();
    //array that stores each of the counts
    int[] give = new int[temp.length+1];

    //copies temp values in ret
    give[0] = totalCount();
    for(int i = 0; i < temp.length; i++){
      give[i+1] = temp[i];
    }
    return give;


  }
  //This will print the current counts for the amino acids
  public static void printAminoAcidCounts(int[] array){

    for(int i = 0; i < array.length; i++){
      System.out.print(array[i]);
    }
    System.out.println();
  }

  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){
    if(next == null){
      return true;
    }

    if(aminoAcid > next.aminoAcid)
      return false;

    return next.isSorted();

}


  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence){
    AminoAcidLL list = new AminoAcidLL(inSequence.substring(0,3));
    inSequence = inSequence.substring(3);
    while(inSequence.length() > 2) {
      list.addCodon(inSequence.substring(0,3));
      inSequence = inSequence.substring(3);
    }
    return list;
  }



  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList){
    AminoAcidLL count = inList;
    char temp;
    while(count.next != null) {
      if(count.aminoAcid > count.next.aminoAcid) {
        temp = count.aminoAcid;
        count.aminoAcid = count.next.aminoAcid;
        count.next.aminoAcid = temp;
        count = count.next;

        if(count.next == null) {
          if(inList.isSorted())
            break;
          count = inList;
        }
      }
      else {
        count = count.next;
        if(count.next == null) {
          if(inList.isSorted())
            break;
          count = inList;
        }
      }
    }
    return inList;
  }
}