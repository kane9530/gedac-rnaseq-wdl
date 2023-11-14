version 1.0

workflow countUniqueItems{
    input {
         Array[String] input_array
    } 

    call countArrayUniqueItems{
        input:
            input_array = input_array
    }

    output{
        Int numberUniqueElements = countArrayUniqueItems.numberUniqueElements
    }
}


task countArrayUniqueItems{
    meta {
        description: "Counts the number of unique elements within an array."
        author: "Kane Toh"
        email: "kanetoh@nus.edu.sg" 
  }

  input {
    Array[String] input_array
  }

  runtime {
    docker: "ubuntu:latest"
  }

  command <<<
    set -eox pipefail
    # WDL requires us to parse in the array with a separator. Recreate the array by
    # wrapping output in ()
    inputArray=~{sep=',' input_array}
    # Command to count number of unique elements and assign to variable
    numberUniqueElements="$(echo "${inputArray}" | tr "," "\n" | sort -u | wc -l)"
    # Print to stdout()
    echo "$numberUniqueElements"
  >>>

  output {
    Int numberUniqueElements = read_int(stdout())
  }

}