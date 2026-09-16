'''Read MCMC inversion binary sample files.'''

from pathlib import Path

import numpy as np
import pandas as pd

# *----------------------------------------------------------------------------* #

SAMPLE_MAGIC = b'SMPL'
NUM_LINES = 20
NUM_REGIONS = 20
MAX_MODEL_PARAMS = 4+3*NUM_LINES+2*NUM_REGIONS
SAMPLE_HEADER_INTS = MAX_MODEL_PARAMS+23
SAMPLE_META_OFFSET = 7+MAX_MODEL_PARAMS
INT_SIZE = 4
FLOAT_SIZE = 4

# *----------------------------------------------------------------------------* #

def parameter_name(model_index, nlines, nregions, magnetic_mode):
    '''
    Purpose
    ----------
    Return the name of a model parameter stored in the sample file.

    Parameters
    ----------
    model_index : int
        Parameter index in the complete inversion model.
    nlines : int
        Number of spectral lines in the model.
    nregions : int
        Number of wavelength regions in the model.
    magnetic_mode : int
        Magnetic representation: 0 for Cartesian and 1 for spherical.

    Returns
    -------
    str
        Human-readable parameter name.

    Notes
    -----
    Unknown indices are returned as ``parameterN``.
    '''

    if model_index < 0:
        return f'parameter{model_index}'

    if model_index < 4:
        magnetic = (
            ['Bz', 'Bx', 'By']
            if magnetic_mode == 0
            else ['Bmod', 'ThetaB', 'PhiB']
        )
        return magnetic[model_index] if model_index < 3 else 'Vlos'

    line_end = 4 + 3 * nlines
    if model_index < line_end:
        line, kind = divmod(model_index - 4, 3)
        return f"line{line}_{['Dopp', 'Damp', 'Eta'][kind]}"

    region_offset = model_index - line_end
    if 0 <= region_offset < 2 * nregions:
        region, kind = divmod(region_offset, 2)
        return f"region{region}_{['Continuum', 'Beta'][kind]}"

    return f'parameter{model_index}'

# *----------------------------------------------------------------------------* #

class SampleFile:
    '''
    Purpose
    ----------
    Read a fixed-record ``Samples.bin`` file written by StokesInversion.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the binary sample file.

    Notes
    -----
    The file uses native-endian 32-bit integers and floating-point values.
    Sample arrays are accessed with memory mapping when a pixel is read.
    '''

    def __init__(self, path):

      self.path = Path(path)
      if not self.path.is_file():
        raise FileNotFoundError(self.path.resolve())

      self.int_dtype = np.dtype('=i4')
      self.float_dtype = np.dtype('=f4')
      self.header = np.fromfile(
          self.path, dtype=self.int_dtype, count=SAMPLE_HEADER_INTS)
      
      if self.header.size != SAMPLE_HEADER_INTS:
        raise ValueError('Truncated sample header')
      
      if self.header[:1].tobytes() != SAMPLE_MAGIC:
        raise ValueError(
            f'Invalid sample magic: {self.header[:1].tobytes()!r}')

      self.nx, self.ny = map(int, self.header[1:3])
      self.nparams = int(self.header[3])
      self.max_samples = int(self.header[4])
      self.float_bits = int(self.header[5])
      self.nmodel = int(self.header[6])

      if self.nx < 1 or self.ny < 1 or self.max_samples < 1:
        raise ValueError('Invalid sample dimensions')
      
      if not 1 <= self.nmodel <= MAX_MODEL_PARAMS:
        raise ValueError('Invalid model-parameter count')
      
      if self.float_bits != 32 or not 0 <= self.nparams <= self.nmodel:
        raise ValueError('Invalid parameter count or sample precision')

      self.model_indices = self.header[7:7+self.nparams].astype(int)

      if np.any((self.model_indices < 0) 
          | (self.model_indices >= self.nmodel)):
        raise ValueError('Invalid model-parameter index')

      meta = self.header[SAMPLE_META_OFFSET : SAMPLE_META_OFFSET + 9]
      self.x_begin, self.x_end, self.y_begin, self.y_end = map(int, meta[:4])

      hash_words = np.asarray(meta[4:6], dtype=self.int_dtype)

      self.config_hash = int(hash_words.view(np.dtype('=u8'))[0])

      self.nlines, self.nregions, self.magnetic_mode = map(int, meta[6:9])

      if not 0 <= self.nlines <= NUM_LINES:
        raise ValueError('Invalid spectral-line count')
      
      if not 0 <= self.nregions <= NUM_REGIONS:
        raise ValueError('Invalid wavelength-region count')
      
      if self.magnetic_mode not in (0, 1):
        raise ValueError('Invalid magnetic mode')
      
      self.names = [parameter_name(index, self.nlines, self.nregions, 
          self.magnetic_mode) for index in self.model_indices]

      self.header_bytes = SAMPLE_HEADER_INTS*INT_SIZE
      self.record_bytes = (2*INT_SIZE+self.nparams*self.max_samples
          *FLOAT_SIZE)
      
      expected_size = self.header_bytes+self.nx*self.ny*self.record_bytes
      actual_size = self.path.stat().st_size

      if actual_size != expected_size:
        raise ValueError(
            f'File-size mismatch: expected {expected_size:,}, '
            f'found {actual_size:,} bytes')

    def _pixel_index(self, x, y):
      '''
      Purpose
      ----------
      Convert an image coordinate to a zero-based file-record index.

      Parameters
      ----------
      x, y : int
          Pixel coordinates in the original observation image.

      Returns
      -------
      int
          Linear pixel-record index.

      Notes
      -----
      Coordinates use the inclusive inversion bounds stored in the header.
      '''

      if not (self.x_begin <= x <= self.x_end 
          and self.y_begin <= y <= self.y_end):
        raise IndexError(f'Pixel {(x, y)} is outside the inversion box')

      ix = x-self.x_begin
      iy = y-self.y_begin
      if ix >= self.nx or iy >= self.ny:
        raise IndexError('Pixel metadata is inconsistent with file dimensions')
      
      return iy*self.nx+ix

    def record_info(self, x, y):
      '''
      Purpose
      ----------
      Read the completion state and sample count for one pixel.

      Parameters
      ----------
      x, y : int
          Pixel coordinates in the original observation image.

      Returns
      -------
      tuple of (bool, int)
          Completion flag and number of valid samples.

      Notes
      -----
      This method reads only the two-integer record header.
      '''

      offset = self.header_bytes+self._pixel_index(x, y)*self.record_bytes
      record = np.fromfile(
          self.path, dtype=self.int_dtype, count=2, offset=offset)
      
      if record.size != 2:
        raise ValueError(f'Truncated sample record for pixel {(x, y)}')
      
      complete, nsamples = record
      return bool(complete), int(nsamples)

    def read_pixel(self, x, y):
      '''
      Purpose
      ----------
      Read all stored samples for one completed pixel.

      Parameters
      ----------
      x, y : int
          Pixel coordinates in the original observation image.

      Returns
      -------
      pandas.DataFrame
          Samples arranged as rows with one column per stored parameter.

      Notes
      -----
      The on-disk parameter-major array is memory mapped and then converted
      to a standalone DataFrame in sample-major order.
      '''

      pixel_index = self._pixel_index(x, y)
      complete, nsamples = self.record_info(x, y)

      if not complete:
        raise ValueError(f'Pixel {(x, y)} is incomplete')
      if not 0 < nsamples <= self.max_samples:
        raise ValueError(f'Invalid sample count {nsamples} for pixel {(x, y)}')
      
      offset = self.header_bytes + pixel_index * self.record_bytes + 2 * INT_SIZE
      raw = np.memmap(self.path, dtype=self.float_dtype, mode='r', 
          offset=offset, shape=(self.nparams, self.max_samples))
      return pd.DataFrame(
          np.asarray(raw[:, :nsamples]).T.copy(), columns=self.names)


    def __repr__(self):
      '''
      Purpose
      ----------
      Return a concise description of the sample file.

      Parameters
      ----------
      None

      Returns
      -------
      str
          Summary of the path, dimensions, parameters, and model metadata.

      Notes
      -----
      The configuration hash is displayed as a 16-digit hexadecimal value.
      '''

      mode = 'CARTESIAN' if self.magnetic_mode == 0 else 'SPHERICAL'
      return (
        f'SampleFile(path={str(self.path)!r}, image={self.nx}x{self.ny}, '
        f'box=({self.x_begin}:{self.x_end}, '
        f'{self.y_begin}:{self.y_end}), nparams={self.nparams}, '
        f'max_samples={self.max_samples}, nlines={self.nlines}, '
        f'nregions={self.nregions}, mode={mode}, '
        f'config_hash=0x{self.config_hash:016x})')

# *----------------------------------------------------------------------------* #
