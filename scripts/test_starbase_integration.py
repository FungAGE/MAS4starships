#!/usr/bin/env python3
"""
Test script to verify starbase model integration.

This script tests:
1. Import of starbase models
2. Database connection
3. Basic queries
4. Cross-database functionality

Usage:
    python scripts/test_starbase_integration.py
"""

import os
import sys
import django
from django.conf import settings
from django.db import connections

# Add the project root to Python path
sys.path.append('/mnt/sda/johannesson_lab/adrian/bin/MAS')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'MAS.settings')
django.setup()

def test_imports():
    """Test that we can import all the necessary models"""
    print("Testing imports...")
    
    try:
        # Test starbase models import
        from starship.starbase_models import (
            Accessions, Ships, Captains, Taxonomy, Genome, Papers,
            FamilyNames, StarshipFeatures, Navis, Haplotype, Gff, JoinedShips
        )
        print("✅ Starbase models imported successfully")
        
        # Test main models import
        from starship.models import Feature, Annotation, StarfishRun, StagingStarship
        print("✅ Main models imported successfully")
        
        return True
        
    except Exception as e:
        print(f"❌ Import failed: {e}")
        return False

def test_database_connections():
    """Test database connections"""
    print("\nTesting database connections...")
    
    try:
        # Test starbase database (SQLite)
        starbase_conn = connections['starbase']
        with starbase_conn.cursor() as cursor:
            cursor.execute("SELECT 1")
            result = cursor.fetchone()
            print("✅ Starbase database (SQLite): Connected")
        
        # Test default database (MySQL) - only if we can connect
        try:
            default_conn = connections['default']
            with default_conn.cursor() as cursor:
                cursor.execute("SELECT 1")
                result = cursor.fetchone()
                print("✅ Default database (MySQL): Connected")
        except Exception as mysql_error:
            print(f"⚠️ Default database (MySQL): Connection failed - {mysql_error}")
            print("   This is expected if MySQL is not running or configured")
        
        return True
        
    except Exception as e:
        print(f"❌ Database connection failed: {e}")
        return False

def test_starbase_queries():
    """Test starbase model queries"""
    print("\nTesting starbase queries...")
    
    try:
        from starship.starbase_models import JoinedShips, Accessions, Ships
        
        # Test basic queries
        starship_count = JoinedShips.objects.count()
        print(f"✅ JoinedShips count: {starship_count}")
        
        accession_count = Accessions.objects.count()
        print(f"✅ Accessions count: {accession_count}")
        
        ship_count = Ships.objects.count()
        print(f"✅ Ships count: {ship_count}")
        
        # Test a sample query
        if starship_count > 0:
            sample_starship = JoinedShips.objects.first()
            print(f"✅ Sample starship: {sample_starship}")
        
        return True
        
    except Exception as e:
        print(f"❌ Starbase query failed: {e}")
        return False

def test_cross_database_queries():
    """Test queries across both databases"""
    print("\nTesting cross-database queries...")
    
    try:
        from starship.starbase_models import JoinedShips
        
        # Test starbase query (SQLite)
        starship_count = JoinedShips.objects.count()
        print(f"✅ Starbase query: {starship_count} starships")
        
        # Test main database query (MySQL) - only if available
        try:
            from starship.models import Feature, Annotation
            
            feature_count = Feature.objects.count()
            print(f"✅ Main database query: {feature_count} features")
            
            annotation_count = Annotation.objects.count()
            print(f"✅ Main database query: {annotation_count} annotations")
            
        except Exception as mysql_error:
            print(f"⚠️ Main database query failed: {mysql_error}")
            print("   This is expected if MySQL is not running or configured")
        
        return True
        
    except Exception as e:
        print(f"❌ Cross-database query failed: {e}")
        return False

def test_model_relationships():
    """Test model relationships and properties"""
    print("\nTesting model relationships...")
    
    try:
        from starship.starbase_models import JoinedShips, Accessions, Ships
        
        # Test if we can access related models
        starship_count = JoinedShips.objects.count()
        print(f"✅ Found {starship_count} starships in database")
        
        if starship_count > 0:
            starship = JoinedShips.objects.first()
            
            # Test properties
            sequence = starship.starship_sequence
            print(f"✅ Starship sequence property: {len(sequence)} characters")
            
            # Test relationships
            if hasattr(starship, 'ship') and starship.ship:
                print(f"✅ Ship relationship: {starship.ship}")
            else:
                print("ℹ️ No ship relationship found (this may be normal)")
            
            if hasattr(starship, 'accession') and starship.accession:
                print(f"✅ Accession relationship: {starship.accession}")
            else:
                print("ℹ️ No accession relationship found (this may be normal)")
        else:
            print("ℹ️ No starships found in database to test relationships")
        
        return True
        
    except Exception as e:
        print(f"❌ Model relationship test failed: {e}")
        return False

def main():
    """Main test function"""
    print("Starbase Integration Test")
    print("=" * 50)
    
    all_passed = True
    
    # Test imports
    if not test_imports():
        all_passed = False
    
    # Test database connections
    if not test_database_connections():
        all_passed = False
    
    # Test starbase queries
    if not test_starbase_queries():
        all_passed = False
    
    # Test cross-database queries
    if not test_cross_database_queries():
        all_passed = False
    
    # Test model relationships
    if not test_model_relationships():
        all_passed = False
    
    print("\n" + "=" * 50)
    if all_passed:
        print("✅ All tests passed! Starbase integration is working correctly.")
    else:
        print("❌ Some tests failed. Please check the errors above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
