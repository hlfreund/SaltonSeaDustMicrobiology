## This function will turn your site x species matrix into a presence/absence table!

counts_to_binary <- function(dataFrame){
  new_m <- matrix(nrow=dim(dataFrame)[1],ncol = dim(dataFrame)[2]) # create new matrix w/ same rows and cols as input dataframe
  ## dim(df)[1] gives you first dimensions (x aka rows), dim(df)[2] gives you second dimensions (y aka columns)

  for( currentRow in 1:nrow(dataFrame)){ # for every row
    for( currentCol in 1:ncol(dataFrame)){ # for every column

      if ( is.na(dataFrame[currentRow, currentCol]) & is.numeric(dataFrame[currentRow, currentCol])){ # if both row and col (specifies each cell) are NA, change val to 0
        new_m[currentRow, currentCol] = 0
        # is.numeric(df[currentRow,currentCol]) is to confirm each cell contains a numeric element
      } else if( is.numeric(dataFrame[currentRow, currentCol]) & dataFrame[currentRow, currentCol] > 0){ # if both row and col (specifies each cell) are > 0, change val to 1
        new_m[currentRow, currentCol] = 1
      } else if ( is.numeric(dataFrame[currentRow, currentCol]) & dataFrame[currentRow, currentCol] == 0){ # if both row and col (specifies each cell) == 0 , change val to 0
        new_m[currentRow, currentCol] = 0
      } else if ( is.character(dataFrame[currentRow, currentCol])){ # if both row and col (specifies each cell) == 0 , change val to 0
        new_m[currentRow, currentCol] = dataFrame[currentRow, currentCol]
      }
    }
  }
  new_df <- as.data.frame(new_m) #turns matrix into dataframe
  names(new_df) <- names(dataFrame) #names rows & cols of new dataframe to be same as row names and col names from input dataframe
  rownames(new_df) <- rownames(dataFrame)
#  new_df2=new_df[,order(ncol(new_df):1)]
  new_df2=new_df[rownames(dataFrame),colnames(dataFrame)]
  return(new_df2) # ensures only output is the new dataframe
}

## Function usage looks like this
# new_df2 <-counts_to_binary(dataFrame)

## Can confirm if new dataframe rows and columns match input dataframe with identical()
## Try: identical (colnames(name of new_df2),colnames(name of dataFrame))
### new_df2 -- whatever you call your new dataframe; dataFrame -- your input dataframe for the function


#new_m <- matrix(nrow=dim(test)[1],ncol = dim(test)[2]) # create new matrix w/ same rows and cols as input dataframe
## dim(df)[1] gives you first dimensions (x aka rows), dim(df)[2] gives you second dimensions (y aka columns)
#
# for( currentRow in 1:nrow(test)){ # for every row
#   for( currentCol in 1:ncol(test)){ # for every column
#
#     if ( is.na(test[currentRow, currentCol]) & is.numeric(test[currentRow, currentCol])){ # if both row and col (specifies each cell) are NA, change val to 0
#       new_m[currentRow, currentCol] = 0
#       # is.numeric(df[currentRow,currentCol]) is to confirm each cell contains a numeric element
#     } else if( is.numeric(test[currentRow, currentCol]) & test[currentRow, currentCol] > 0){ # if both row and col (specifies each cell) are > 0, change val to 1
#       new_m[currentRow, currentCol] = 1
#     } else if ( is.numeric(test[currentRow, currentCol]) & test[currentRow, currentCol] == 0){ # if both row and col (specifies each cell) == 0 , change val to 0
#       new_m[currentRow, currentCol] = 0
#     } else if ( is.character(test[currentRow, currentCol])){ # if both row and col (specifies each cell) == 0 , change val to 0
#       new_m[currentRow, currentCol] = test[currentRow, currentCol]
#     }
#   }
# }
# new_df <- as.data.frame(new_m) #turns matrix into dataframe
# names(new_df) <- names(test) #names rows & cols of new dataframe to be same as row names and col names from input dataframe
# rownames(new_df) <- rownames(test)
# #new_df2=new_df[,order(ncol(new_df):1)]
# new_df2=new_df[rownames(test),colnames(test)]
