import argparse
import os
import subprocess
import concurrent.futures
import xml.etree.ElementTree as ET
import tempfile
import re
import shutil

namespace = {"ns": "http://regis-web.systemsbiology.net/pepXML"}

img_pattern = re.compile(r'show_tmp_pngfile\.pl\?file=(.+?\.png)')
err_pattern = re.compile(r'([Ee]rror .+)')
ver_pattern = re.compile(r'(TPP v[^<]+)')

def start_processing(input_path, mzml_path, parser_path, output_path, img_dir, tsv_path, num_jobs, debug):

    print("Start processing...")

    had_errors = False

    if tsv_path:
        tsv_file = open(tsv_path, 'w')
        tsv_file.write(f"spectrum\trank\told_decimal_ratio\tlight_area\theavy_area\tdecimal_ratio\n")
    else:
        tsv_file = None

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_jobs) as executor:
        # Parse input XML file
        ET.register_namespace("", namespace["ns"])
        tree = ET.parse(input_path)
        root = tree.getroot()
        xpressratio_summary = root.find('.//ns:xpressratio_summary', namespace)
        spectrum_queries = root.findall('.//ns:spectrum_query', namespace)

        futures = []
        for spectrum_query in spectrum_queries:
            spectrum = spectrum_query.get('spectrum')
            # Search for nested xpressratio_result elements
            xpress_results = spectrum_query.findall('.//ns:analysis_result[@analysis="xpress"]/ns:xpressratio_result', namespace)

            # Process each xpressratio_result in parallel
            for i, xpress_result in enumerate(xpress_results):
                future = executor.submit(process_xpress_result, xpress_result, spectrum_query, xpressratio_summary, mzml_path, parser_path, img_dir, debug, i + 1)
                futures.append(future)

        jobs_total = len(futures)
        jobs_done = 0
        for i, future in enumerate(concurrent.futures.as_completed(futures)):
            try:
                spectrum, rank, old_decimal_ratio, light_area, heavy_area, decimal_ratio = future.result()
                jobs_done += 1
                jobs_pct = round(jobs_done/jobs_total*100, 1)
                print(f"{spectrum} - ratio updated from {old_decimal_ratio} to {decimal_ratio} - {jobs_pct}%")
                if tsv_file:
                    tsv_file.write(f"{spectrum}\t{rank}\t{old_decimal_ratio}\t{light_area}\t{heavy_area}\t{decimal_ratio}\n")
            except Exception as e:
                print(f"Error: {e}")
                had_errors = True
                executor.shutdown(wait=False)
                for f in futures:
                    if not f.done():
                        f.cancel()
                if tsv_file:
                    tsv_file.close()
                break

    if not had_errors:
        try:
            tree.write(output_path, encoding='utf-8', xml_declaration=True)
            print("Processing complete.")
        except Exception as e:
            print(f"Error: {e}")

    if tsv_file and not tsv_file.closed:
        tsv_file.close()

def process_xpress_result(xpress_result, spectrum_query, xpressratio_summary, mzml_path, parser_path, img_dir, debug, rank):
    # Extract relevant attributes from xpress_result
    light_firstscan = xpress_result.get('light_firstscan')
    light_lastscan = xpress_result.get('light_lastscan')
    heavy_firstscan = xpress_result.get('heavy_firstscan')
    heavy_lastscan = xpress_result.get('heavy_lastscan')
    light_mass = xpress_result.get('light_mass')
    heavy_mass = xpress_result.get('heavy_mass')
    mass_tol = xpress_result.get('mass_tol')
    old_decimal_ratio = xpress_result.get('decimal_ratio')

    # Extract ppmtol, min_num_isotope_peaks, and xpress_light from xpressratio_summary
    ppmtol = xpressratio_summary.get('ppmtol')
    min_num_isotope_peaks = xpressratio_summary.get('min_num_isotope_peaks')
    xpress_light = xpressratio_summary.get('xpress_light')

    # Extract spectrum attributes from the parent spectrum_query tag
    spectrum = spectrum_query.get('spectrum')
    assumed_charge = spectrum_query.get('assumed_charge')
    index = spectrum_query.get('index')

    if img_dir:
        os.makedirs(img_dir, exist_ok=True)

    light_area = heavy_area = decimal_ratio = None
    # Create a temporary directory to store the reduced XML file
    with tempfile.TemporaryDirectory() as temp_dir:
        reduced_xml_path = os.path.join(temp_dir, "reduced_xml_to_modify.xml")
        mzml_link_path = os.path.join(temp_dir, "spectra.mzML")
        os.symlink(os.path.realpath(mzml_path), mzml_link_path)

        # Create a reduced XML file containing only the necessary elements
        reduced_xml = ET.Element("msms_pipeline_analysis")
        msms_run_summary = ET.SubElement(reduced_xml, "msms_run_summary")
        spectrum_query = ET.SubElement(msms_run_summary, "spectrum_query", spectrum=spectrum, assumed_charge=assumed_charge, index=index)
        search_hit = ET.SubElement(spectrum_query, "search_hit", hit_rank = "1")
        analysis_result = ET.SubElement(search_hit, "analysis_result", analysis="xpress")
        
        # Copy the entire xpress_result element to the new xpressratio_result
        xpressratio_result = ET.SubElement(analysis_result, "xpressratio_result")
        for attrib in xpress_result.attrib:
            xpressratio_result.set(attrib, xpress_result.get(attrib))

        # Write the reduced XML to a file
        reduced_tree = ET.ElementTree(reduced_xml)
        reduced_tree.write(reduced_xml_path)

        # Prepare parameters for XPressPeptideUpdateParser
        query = {
            'LightFirstScan': light_firstscan,
            'LightLastScan': light_lastscan,
            'HeavyFirstScan': heavy_firstscan,
            'HeavyLastScan': heavy_lastscan,
            'XMLFile': mzml_link_path,
            'ChargeState': assumed_charge,
            'LightMass': light_mass,
            'HeavyMass': heavy_mass,
            'MassTol': mass_tol,
            'PpmTol': ppmtol,
            'NumIsotopePeaks': min_num_isotope_peaks,
            'index': index,
            'xmlfile': reduced_xml_path,
            'FileOut': spectrum, # TPP v. 5
            'OutFile': spectrum, # TPP v. 6
            'bXpressLight1': xpress_light
        }
        query['overwrite'] = '1'
        query['no_images'] = '0' if img_dir else '1' # for custom versions supporting this
        query_string = '&'.join(['{}={}'.format(k,v) for k,v in query.items()])

        # Call XPressPeptideUpdateParser
        result = subprocess.run([parser_path], env={'QUERY_STRING': query_string, 'REQUEST_METHOD': 'GET'}, capture_output=True, text=True, check=True)
        # Check if the result contains 'ratio updated'
        if (debug):
            print(result.stdout)
        if 'ratio updated' in result.stdout or 'Ratio updated' in result.stdout:
            # Parse the modified reduced XML file
            modified_tree = ET.parse(reduced_xml_path)
            modified_root = modified_tree.getroot()
            modified_xpress_result = modified_root.find('.//xpressratio_result')

            decimal_ratio = modified_xpress_result.get('decimal_ratio')
            light_area = modified_xpress_result.get('light_area')
            heavy_area = modified_xpress_result.get('heavy_area')

            # Update the original XML file with the modified values
            for attrib in modified_xpress_result.attrib:
                xpress_result.set(attrib, modified_xpress_result.get(attrib))
            xpress_result.set("heavy2light_ratio", f"1:{decimal_ratio}")

            if img_dir:
                light_png_src, heavy_png_src = img_pattern.findall(result.stdout)
                light_png_dest = os.path.join(img_dir, f"{spectrum}_{rank}_0.png")
                heavy_png_dest = os.path.join(img_dir, f"{spectrum}_{rank}_1.png")
                shutil.move(light_png_src, light_png_dest)
                shutil.move(heavy_png_src, heavy_png_dest)
        else:
            errors = err_pattern.findall(result.stdout)
            raise Exception(f"error running XPressPeptideUpdateParser: {errors}")

    return spectrum, rank, old_decimal_ratio, light_area, heavy_area, decimal_ratio

def find_parser_path():

    # Check if XPressPeptideUpdateParser is in the same directory as the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidate = os.path.join(script_dir, "XPressPeptideUpdateParser.cgi")
    if os.path.isfile(candidate):
        return candidate

    # Check if XPressPeptideUpdateParser is available in PATH
    for path in os.environ["PATH"].split(os.pathsep):
        candidate = os.path.join(path.strip('"'), "XPressPeptideUpdateParser.cgi")
        if os.path.isfile(candidate):
            return candidate

    return None

def main():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description='Run XPressPeptideUpdateParser batch')
    parser.add_argument('--input', help='Input pep.xml file')
    parser.add_argument('--mzml', help='mzML file')
    parser.add_argument('--output', help='Output pep.xml file')
    parser.add_argument('--tsv', help='Save ratios in the tsv file')
    parser.add_argument('--img', help='Directory to save the images if wanted')
    parser.add_argument('--parser', help='Path to XPressPeptideUpdateParser.cgi')
    parser.add_argument('--jobs', default=2, type=int, help='Number of parallel jobs [2]')
    parser.add_argument('--debug', default=False, action='store_true', help='Debugging mode')

    # Parse command-line arguments
    args = parser.parse_args()

    if not args.input:
        parser.error('Input pep.xml file not specified')
    if not args.mzml:
        parser.error('mzML file not specified')
    if not args.output:
        parser.error('Output pep.xml file not specified')
    if not args.parser:
        args.parser = find_parser_path()
        if not args.parser:
            parser.error('Path to XPressPeptideUpdateParser.cgi not specified')

    if not os.path.isfile(args.input):
        parser.error('Specified input pep.xml file does not exist')
    if not os.path.isfile(args.mzml):
        parser.error('Specified mzML file does not exist')
    if not os.path.isfile(args.parser):
        parser.error('Specified path to XPressPeptideUpdateParser.cgi does not exist')

    result = subprocess.run([args.parser], env={'REQUEST_METHOD': 'GET'}, capture_output=True, text=True)
    version = ver_pattern.findall(result.stdout)
    if not version:
        parser.error(f'Could not run XPressPeptideUpdateParser.cgi, got:\n{result.stdout}')
    version = version[0].strip(' :').replace('))', ')')
    print(f"Using XPressPeptideUpdateParser.cgi from {version}")

    # Start processing
    start_processing(args.input, args.mzml, args.parser, args.output, args.img, args.tsv, int(args.jobs), args.debug)

if __name__ == "__main__":
    main()
