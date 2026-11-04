# CHIMERA v2.2 benchmark-data request (draft; not sent)

**Status:** Draft for the project maintainer to review and send. No message has been sent by this repository workflow.

## Subject
Request for reproducible input data for the CHIMERA v2.2 incomplete-galaxy-catalogue example

## Message

Hello CHIMERA team,

I am working on a small, reproducibility-focused methodological study of galaxy-catalogue incompleteness and host-galaxy weighting in dark-siren inference. I am trying to reproduce the incomplete-catalogue example using the exact pinned v2.2 notebook, rather than substituting a different mock and describing it as equivalent.

The notebook at commit `fd35f897e36850d865c760c13325c5a46ebe94db` appears to require the following data products, which I have not been able to locate in the public release/repository:

- `O5snr25_100evs_mock1_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5`
- `injections_Ninj_2e7_O5Like_snr20.h5`
- `completeness_masks_5x_pix32False.h5`
- `completeness_masks_5x_pix32False_mock1_mcut10.h5`
- `temp_gw_mask_idxs_mock1_mcut10.h5`
- `pcat_incomplete.h5`

Could you point me to an approved download/archive or repository revision containing the exact files and any accompanying manifest, or advise whether these products are available on request? If redistribution is restricted, information about the access procedure and applicable licence/terms would be very helpful. I can provide a file manifest and SHA-256 checksums in the resulting reproducibility record.

I will not report the example as reproduced unless the inputs can be traced to an authoritative source and the pinned workflow can be executed and audited. I understand if the data are not publicly distributable.

Thank you for maintaining CHIMERA and for any guidance you can provide.

Best regards,  
Vedant

## Before sending

- Confirm the recipient/contact route from the official CHIMERA repository.
- Verify the pinned notebook and filenames again against the current intended benchmark.
- Do not attach unpublished or private material.
- Record the date, recipient, response, access terms, and checksums if the request is sent.
- This draft is not evidence that any maintainer has been contacted.
