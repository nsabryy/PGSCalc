############################### SCORE JACCARD INDEX LOGIC 
        # # Initialize an empty DataFrame with columns for chromosome, position, score count, and Jaccard index.
        # snp_score_count = pd.DataFrame(columns=["chromosome", "position", "score_count", "jaccard_index"])

        # # Collect all unique traits (pgs_id) for the union set.
        # total_scores = len(config["PGScore"]["pgs_ids"])

        # # Iterate over the combined_scores dictionary.
        # for chrom, pos_list in score_package.combined_scores.items():
        #     for pos, scores in pos_list.items():
        #         # Extract traits (pgs_ids) associated with this SNP.
        #         traits = {score['score'] for score in scores}
                
        #         # Calculate the Jaccard index for this SNP.
        #         jaccard_index = len(traits) / len(all_traits) if all_traits else 0

        #         # Append the information to the DataFrame.
        #         snp_score_count = snp_score_count.concat(
        #             {
        #                 "chromosome": chrom,
        #                 "position": pos,
        #                 "score_count": len(scores),
        #                 "jaccard_index": jaccard_index
        #             }, 
        #             ignore_index=True
        #         )

        # # Save the DataFrame to a CSV file.
        # snp_score_count.to_csv("/N/project/compgen/PGSCalc/scoring_results/snp_score_count.csv", index=False)
        # comm.Abort() 

###############################