class DatabaseRouter:
    """
    Database router to handle routing between the main 'mas' database 
    and the 'starbase' database for starbase-related models.
    """
    
    # Models that should use the starbase database
    STARBASE_MODELS = {
        'accessions', 'ships', 'captains', 'taxonomy', 'genome', 
        'papers', 'familynames', 'starshipfeatures', 'navis', 
        'haplotype', 'gff', 'joinedships'
    }
    
    def db_for_read(self, model, **hints):
        """Route read operations to appropriate database."""
        if model._meta.model_name in self.STARBASE_MODELS:
            return 'starbase'
        return 'default'

    def db_for_write(self, model, **hints):
        """Route write operations to appropriate database."""
        if model._meta.model_name in self.STARBASE_MODELS:
            return 'starbase'
        return 'default'

    def allow_relation(self, obj1, obj2, **hints):
        """Allow relations between objects in the same database."""
        db1 = self.db_for_read(obj1.__class__)
        db2 = self.db_for_read(obj2.__class__)
        return db1 == db2

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        """Control which database migrations are applied to."""
        if app_label == 'starship':
            # Starbase models should only be migrated to starbase database
            if model_name in self.STARBASE_MODELS:
                return db == 'starbase'
            # Other starship models (like Annotation, Feature) go to default database
            else:
                return db == 'default'
        # All other apps use default database
        return db == 'default'
